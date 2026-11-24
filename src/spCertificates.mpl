# Return the reduced gens and its indices in the caller's original gens.
local remove_redundant_constants := proc(gens)
local i, indices := select(i -> not (type(gens[i], rational) and gens[i] >= 0),
                        [seq(i, i=1..nops(gens))]);
    return map(i -> gens[i], indices), indices;
end proc;

# max_search_time is a CPU-time budget enforced by one timelimit around the
# whole pipeline. Stages need no deadline of their own; their iteration caps
# bound the work, and only the monogenic per-query limit affects scheduling.
spCertificates := proc(f, gens0, x,
                       {max_k::nonnegint := 100, max_bisection::nonnegint := 200,
                        query_time_limit::positive := 1/10, max_query_retries::nonnegint := 4,
                        max_search_time::positive := 300,
                        use_univariate_heuristics::truefalse := true,
                        max_heuristic_trials::nonnegint := 64,
                        max_scalar_updates::nonnegint := 32,
                        max_averkov_N::nonnegint := 256})
# The deadline is a local passed down explicitly: module state would go
# stale when a caller's own timelimit skips this routine's finally block.
local stage := "preprocessing", search_deadline, pipeline;
$ifdef LOG_TIME
local entry_stack_level;
    INIT_START_LOG_TIME("spCertificates")
    entry_stack_level := stack_level;
$endif

    pipeline := proc()
    local gens, original_count, indices, remaining;
    local g, f2, monogenic, certificates, i;
$ifdef ENABLE_QUICK_CERT
    local quickcheck, quickcertificates;
$endif
        # Expand before classifying constants or extracting polynomial coefficients.
        if not type(gens0,list) then error "gens must be a list of rational univariate polynomials"; end if;
        original_count := nops(gens0);
        gens := map(expand,gens0);
        if not type(f,polynom(rational,x)) or
        not andmap(g_i->type(g_i,polynom(rational,x)),gens) then
            error "f and gens must be univariate polynomials with rational coefficients";
        end if;
        gens,indices := remove_redundant_constants(gens);
        if univariate_nonnegative_on_R(f,x) then
            return restore_generator_multipliers([f,op(map(g_i->0,gens))],indices,original_count);
        end if;
        if gens = [] then
            error "no nonredundant generators remain and f is not nonnegative on the real line";
        end if;
        for i from 1 to nops(gens) do
            if type(gens[i], rational) and gens[i] < 0 then
                certificates := map(gen -> 0, gens);
                certificates[i] := ((f - 1)/2)^2/(-gens[i]);
                return restore_generator_multipliers([((f + 1)/2)^2,op(certificates)],indices,original_count);
            end if;
        end do;

$ifdef ENABLE_QUICK_CERT
        stage := "quick certificates";
        quickcheck, quickcertificates := quickCertificates(f,gens,x);
        if quickcheck then
            DEBUG(__FILE__, __LINE__, lprint(">> Quick certificates found", op(quickcertificates)));
            return restore_generator_multipliers(quickcertificates,indices,original_count);
        end if;
$endif

        stage := "bound construction";
        try
            g := bound_poly(gens, x);
        catch "no compactness witness found":
            error "bound_poly could not witness compactness of gens: expected an even-degree generator with negative leading coefficient or two odd-degree generators with opposite-sign leading coefficients";
        end try;
        DEBUG(__FILE__, __LINE__, lprint(">> Done with bound_poly"));
        DEBUG(__FILE__, __LINE__, lprint(">> bound_poly g", g));

        DEBUG(__FILE__, __LINE__, lprint(">> Polynomial for compact-bound reduction", f));
        DEBUG(__FILE__, __LINE__, lprint(">> Compact bound", g[1]));
        stage := "univariate reduction";
        certificates := FAIL;
        if use_univariate_heuristics then
            certificates := univariate_reduction(f,gens,x,g[1],
                            ':-max_trials'=max_heuristic_trials, ':-max_updates'=max_scalar_updates);
        end if;
        if certificates=FAIL then
            stage := "Averkov";
            certificates := averkov_lemma_7(x,f,gens,g[1],N_GUESS_AVKL,':-max_N'=max_averkov_N);
        end if;
        DEBUG(__FILE__, __LINE__, lprint(">> Done with compact-bound reduction"));

        DEBUG(__FILE__, __LINE__, lprint(">> Reduction multipliers", certificates));
        f2 := f - dot_product(gens, certificates);
$ifdef ENABLE_AVERKOV_CHECK
        lprint(">> 2. Checking correctness of averkov_lemma_7",
               univariate_sign_set([g[1],f2],["nonnegative","nonpositive"],x));
$endif
        DEBUG(__FILE__, __LINE__, lprint(">> f2", f2));
        stage := "reduction residual";
        if univariate_nonnegative_on_R(f2,x) then
            return restore_generator_multipliers([f2,op(certificates)],indices,original_count);
        end if;

        # Global keyword names avoid substitution of the caller's parameter values.
        stage := "monogenic search";
        remaining := search_deadline-time();
        if remaining <= 0 then error "time expired"; end if;
        monogenic := monogenicCertificates(f2, g[1], x,
                                           ':-max_k'=max_k, ':-max_bisection'=max_bisection,
                                           ':-query_time_limit'=query_time_limit, ':-max_query_retries'=max_query_retries,
                                           ':-max_search_time'=remaining,
                                           ':-search_deadline'=search_deadline);
        DEBUG(__FILE__, __LINE__, lprint(">> Done with monogenicCertificates"));

        certificates[g[4]] := certificates[g[4]] + monogenic[2]*g[2];
        if g[3] <> 0 then
            certificates[g[5]] := certificates[g[5]] + monogenic[2]*g[3];
        end if;

        certificates := [monogenic[1], op(certificates)];
$ifdef ENABLE_VERIFICATION
        lprint(">> This should be zero", expand(f - dot_product([1, op(gens)], certificates)));
$endif
        DEBUG(__FILE__, __LINE__, lprint(">> Input polynomial:", f));
        DEBUG(__FILE__, __LINE__, lprint(">> Generators:", 1, op(gens)));
        DEBUG(__FILE__, __LINE__, lprint(">> Certificates:", op(certificates)));
        return restore_generator_multipliers(certificates,indices,original_count);
    end proc;

    try
        search_deadline := time()+max_search_time;
        return timelimit(max_search_time, pipeline());
    catch "time expired":
$ifdef LOG_TIME
        # Interrupted routines below the timelimit never ran END_LOG_TIME.
        log_unwind(entry_stack_level);
$endif
        error "certificate search exhausted max_search_time during %1",stage;
$ifdef LOG_TIME
    finally
        END_LOG_TIME("spCertificates");
$endif
    end try;
end proc;
