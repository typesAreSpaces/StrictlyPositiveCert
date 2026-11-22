# Return the reduced gens and its indices in the caller's original gens.
local remove_redundant_constants := proc(gens)
local i, indices := select(i -> not (type(gens[i], rational) and gens[i] >= 0),
                        [seq(i, i=1..nops(gens))]);
    return map(i -> gens[i], indices), indices;
end proc;

spCertificates := proc(f, gens, x,
                       {max_k::nonnegint := 100, max_bisection::nonnegint := 200,
                        query_time_limit::positive := 1/10, max_query_retries::nonnegint := 4,
                        max_search_time::positive := 300})
$ifdef LOG_TIME
    INIT_START_LOG_TIME("spCertificates",0)
$endif
    if univariate_nonnegative_on_R(f,x) then
        return [f, op(map(0, gens))];
    end if;
local g, f2, monogenic, certificates, reduced_gens, indices, restored;
local i;
$ifdef ENABLE_QUICK_CERT
local quickcheck, quickcertificates;
$endif
    reduced_gens, indices := remove_redundant_constants(gens);
    if nops(reduced_gens) <> nops(gens) then
        certificates := spCertificates(f, reduced_gens, x,
                                      ':-max_k'=max_k, ':-max_bisection'=max_bisection,
                                      ':-query_time_limit'=query_time_limit, ':-max_query_retries'=max_query_retries,
                                      ':-max_search_time'=max_search_time);
        restored := [certificates[1], op(map(gen -> 0, gens))];
        for i from 1 to nops(indices) do
            restored[indices[i]+1] := certificates[i+1];
        end do;
$ifdef LOG_TIME
        END_LOG_TIME("spCertificates",0);
$endif
        return restored;
    end if;
    if gens = [] then
        error "no nonredundant generators remain and f is not nonnegative on the real line";
    end if;
    for i from 1 to nops(gens) do
        if type(gens[i], rational) and gens[i] < 0 then
            certificates := map(gen -> 0, gens);
            certificates[i] := ((f - 1)/2)^2/(-gens[i]);
            return [((f + 1)/2)^2, op(certificates)];
        end if;
    end do;

$ifdef ENABLE_QUICK_CERT
    quickcheck, quickcertificates := quickCertificates(f, gens, x);
    if quickcheck then
        DEBUG(__FILE__, __LINE__, lprint(">> Quick certificates found", op(quickcertificates)));
        return quickcertificates;
    end if;
$endif

    g := bound_poly(gens, x);
    DEBUG(__FILE__, __LINE__, lprint(">> Done with bound_poly"));
    DEBUG(__FILE__, __LINE__, lprint(">> bound_poly g", g));

    DEBUG(__FILE__, __LINE__, lprint(">> Poly f for averkov_lemma_7", f));
    DEBUG(__FILE__, __LINE__, lprint(">> Poly g[1] for averkov_lemma_7", g[1]));
    certificates := averkov_lemma_7(x, f, gens, g[1], N_GUESS_AVKL);
    DEBUG(__FILE__, __LINE__, lprint(">> Done with averkov_lemma_7"));

    DEBUG(__FILE__, __LINE__, lprint(">> Averkov multipliers", certificates));
    f2 := f - dot_product(gens, certificates);
$ifdef ENABLE_AVERKOV_CHECK
    lprint(">> 2. Checking correctness of averkov_lemma_7",
           univariate_sign_set([g[1],f2],["nonnegative","nonpositive"],x));
$endif
    DEBUG(__FILE__, __LINE__, lprint(">> f2", f2));

    # Global keyword names avoid substitution of the caller's parameter values.
    monogenic := monogenicCertificates(f2, g[1], x,
                                       ':-max_k'=max_k, ':-max_bisection'=max_bisection,
                                       ':-query_time_limit'=query_time_limit, ':-max_query_retries'=max_query_retries,
                                       ':-max_search_time'=max_search_time);
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
$ifdef LOG_TIME
    END_LOG_TIME("spCertificates",0);
$endif
    return certificates;
end proc;
