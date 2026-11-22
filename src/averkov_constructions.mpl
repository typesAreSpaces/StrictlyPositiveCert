# We assume:
# 1. SemiAlgebraic(B_poly) is compact
# Return: list of sums of squares multipliers l
# such that f - dot_product(l, basis) > 0 over SemiAlgebraic(B_poly)
local averkov_lemma_7 := proc(x, f, basis, B_poly, N_guess)
$ifdef LOG_TIME
    INIT_START_LOG_TIME("averkov_lemma_7",0)
$endif
local gamma, interval, lowerbound, upperbound;
local eps;
local N, g, term;
local semialgebraic_of_B;
local M, mu, m, i;
local pos_coeff;
$ifdef ENABLE_DEBUGGING
local T;
$endif
local semialgebraic_for_mu;
local _exp1, _exp2, _exp3;
local pos_coeff1, pos_coeff2, N1, N2;
$ifdef ENABLE_BINARY_SEARCH_AVKL
local N_top, N_bot, N_cur, init_sos, predicate_bin_search;
$endif
    DEBUG(__FILE__, __LINE__, lprint(">> Start @averkov_lemma_7"));
    DEBUG(__FILE__, __LINE__, lprint(">> f", f));
    DEBUG(__FILE__, __LINE__, lprint(">> basis", basis));
    DEBUG(__FILE__, __LINE__, lprint(">> B_poly", B_poly));

$ifdef LOG_TIME
    START_LOG_TIME("averkov_lemma_7::Minimization_f",1);
$endif
    semialgebraic_of_B := univariate_sign_set([B_poly],["nonnegative"],x);
    DEBUG(__FILE__, __LINE__, lprint(">> semialgebraic_of_B", semialgebraic_of_B));

    if evalb(semialgebraic_of_B = []) then
        return map(g_i -> 0, basis);
    end if;

    M := -min(
        map(proc(bound)
                interval := bound_info(x, bound, 0);
                # TOCHECK
                # This might introduce a bug if `lowerbound > upperbound`
                # happens to be true for some reason
                lowerbound := convert(evalf(interval[1]), rational);
                upperbound := convert(evalf(interval[2]), rational);
                simplify(minimize(f, x = lowerbound .. upperbound))
            end proc,
            semialgebraic_of_B)
             );
$ifdef LOG_TIME
    END_LOG_TIME("averkov_lemma_7::Minimization_f",1);
$endif
    DEBUG(__FILE__, __LINE__, lprint(">> M", evalf(M)));
    if evalf(M < 0) then
        DEBUG(__FILE__, __LINE__, lprint(">> Done because f is strictly positive over SemiAlgebraic(B_poly)"));
$ifdef LOG_TIME
        END_LOG_TIME("averkov_lemma_7",0)
$endif
        return map(g_i -> 0, basis);
    end if;

    # Keep the logarithms in the exponent heuristic defined.
    M := max(M, 1);

    m := numelems(basis);

    #
    # Find gamma
    #
$ifdef LOG_TIME
    START_LOG_TIME("averkov_lemma_7::compute_gamma",2);
$endif
    gamma := 1/2*max(
        map(proc(g_i)
                map(
                    proc(bound)
                        interval := bound_info(x, bound, 0);
                        # TOCHECK
                        # This might introduce a bug if `lowerbound > upperbound`
                        # happens to be true for some reason
                        lowerbound := convert(evalf(interval[1]), rational);
                        upperbound := convert(evalf(interval[2]), rational);
                        simplify(maximize(g_i, x = lowerbound .. upperbound))
                    end proc,
                    semialgebraic_of_B)
            end proc,
            basis));
$ifdef LOG_TIME
    END_LOG_TIME("averkov_lemma_7::compute_gamma",2);
$endif
    # We just need a bound, it doesn't need to be
    # the tightest bound [to discuss later]
    gamma := max(ceil(evalf(gamma)), 1);
    DEBUG(__FILE__, __LINE__, lprint(">> gamma", gamma));

    #
    # Find exponent eps
    #
$ifdef LOG_TIME
    START_LOG_TIME("averkov_lemma_7::compute_eps",3);
$endif
$ifdef ENABLE_DEBUGGING
    T := univariate_sign_set([B_poly,f],["nonnegative","negative"],x);
    DEBUG(__FILE__, __LINE__, lprint(">> T", evalf(T)));
$endif

    eps := 1/2*convert(evalf(findEps(x, basis, f)), rational);
    DEBUG(__FILE__, __LINE__, lprint(">> eps", evalf(eps)));
$ifdef LOG_TIME
    END_LOG_TIME("averkov_lemma_7::compute_eps",3);
$endif

    #
    # Find mu
    #
$ifdef LOG_TIME
    START_LOG_TIME("averkov_lemma_7::compute_mu",4);
$endif
    semialgebraic_for_mu := univariate_sign_set([B_poly,op(map(g_i -> g_i+EPS_FACTOR*eps,basis))],
                                                [seq("nonnegative",i=1..nops(basis)+1)],x);
    DEBUG(__FILE__, __LINE__, lprint(">> semialgebraic_for_mu", semialgebraic_for_mu));
    mu := min(
        map(proc(bound)
                interval := bound_info(x, bound, 0);
                # TOCHECK
                # This might introduce a bug if `lowerbound > upperbound`
                # happens to be true for some reason
                lowerbound := convert(evalf(interval[1]), rational);
                upperbound := convert(evalf(interval[2]), rational);
                simplify(minimize(f, x = lowerbound .. upperbound))
            end proc,
            semialgebraic_for_mu)
             );

    mu := convert(evalf(mu), rational);
    DEBUG(__FILE__, __LINE__, lprint(">> mu", mu));
$ifdef LOG_TIME
    END_LOG_TIME("averkov_lemma_7::compute_mu",4);
$endif

$ifdef LOG_TIME
    START_LOG_TIME("averkov_lemma_7::compute_N_heuristic",5);
$endif
    #
    # Find N
    #
    _exp1 := (log(2*m*gamma) - log(alpha*mu))/(log(gamma + eps) - log(gamma));
    _exp2 := (log(2*m*gamma) - log(alpha*M))/(log(gamma + eps) - log(gamma));
    _exp3 := (log(alpha*M) - log(2*eps))/(log(gamma + 2*eps) - log(gamma + eps));
    DEBUG(__FILE__, __LINE__, lprint(">> _exp1", _exp1));
    DEBUG(__FILE__, __LINE__, lprint(">> _exp2", _exp2));
    DEBUG(__FILE__, __LINE__, lprint(">> _exp3", _exp3));

    # 1.
    if evalf(subs(alpha=1, simplify(_exp1) <= simplify(_exp2) and simplify(_exp2) <= simplify(_exp3))) then
        DEBUG(__FILE__, __LINE__, lprint(">> Case 1 @ averkov_lemma_7"));
        pos_coeff1 := convert(
            evalf(
                solve(
                    simplify(_exp2) = simplify(_exp3), alpha, 'maxsols'=1)
                 ),
            rational);
        pos_coeff2 := convert(
            evalf(
                solve(
                    simplify(_exp1) = simplify(_exp3), alpha, 'maxsols'=1)
                 ),
            rational);
        N1 := ceil(1/2*subs(alpha=pos_coeff1, simplify(_exp3)));
        N2 := max(ceil(1/2*subs(alpha=pos_coeff2, simplify(_exp2))), ceil(1/2*subs(alpha=pos_coeff2, simplify(_exp3))));
        if evalf(N1 <= N2) then
            N := N1;
            pos_coeff := pos_coeff1;
        else
            N := N2;
            pos_coeff := pos_coeff2;
        end if;
    end if;

    # 2.
    if evalf(subs(alpha=1, simplify(_exp1) <= simplify(_exp3) and simplify(_exp3) <= simplify(_exp2))) then
        DEBUG(__FILE__, __LINE__, lprint(">> Case 2 @ averkov_lemma_7"));
        pos_coeff := convert(
            evalf(
                solve(
                    simplify(_exp3) = simplify(_exp2), alpha, 'maxsols'=1)
                 ),
            rational);
        N := ceil(1/2*subs(alpha=pos_coeff, simplify(_exp2)));
    end if;

    # 3.
    if evalf(subs(alpha=1, simplify(_exp2) <= simplify(_exp1) and simplify(_exp1) <= simplify(_exp3))) then
        DEBUG(__FILE__, __LINE__, lprint(">> Case 3 @ averkov_lemma_7"));
        pos_coeff1 := convert(
            evalf(
                solve(
                    simplify(_exp1) = simplify(_exp3), alpha, 'maxsols'=1)
                 ),
            rational);
        pos_coeff2 := convert(
            evalf(
                solve(
                    simplify(_exp2) = simplify(_exp3), alpha, 'maxsols'=1)
                 ),
            rational);
        N1 := ceil(evalf(1/2*subs(alpha=pos_coeff1, simplify(_exp3))));
        N2 := max(ceil(evalf(1/2*subs(alpha=pos_coeff2, simplify(_exp1)))), ceil(evalf(1/2*subs(alpha=pos_coeff2, simplify(_exp3)))));
        if N1 <= N2 then
            N := N1;
            pos_coeff := pos_coeff1;
        else
            N := N2;
            pos_coeff := pos_coeff2;
        end if;
    end if;

    # 4.
    if evalf(subs(alpha=1,simplify(_exp2) <= simplify(_exp3) and simplify(_exp3) <= simplify(_exp1))) then
        DEBUG(__FILE__, __LINE__, lprint(">> Case 4 @ averkov_lemma_7"));
        pos_coeff := convert(
            evalf(
                solve(
                    simplify(_exp3) = simplify(_exp1), alpha, 'maxsols'=1)
                 ),
            rational);
        N := ceil(1/2*subs(alpha=pos_coeff, simplify(_exp1)));
    end if;

    # 5.
    if evalf(subs(alpha=1,simplify(_exp3) <= simplify(_exp1) and simplify(_exp1) <= simplify(_exp2))) then
        DEBUG(__FILE__, __LINE__, lprint(">> Case 5 @ averkov_lemma_7"));
        pos_coeff := convert(
            evalf(
                solve(
                    simplify(_exp3) = simplify(_exp2), alpha, 'maxsols'=1)
                 ),
            rational);
        N := ceil(1/2*subs(alpha=pos_coeff, simplify(_exp2)));
    end if;

    # 6.
    if evalf(subs(alpha=1,simplify(_exp3) <= simplify(_exp2) and simplify(_exp2) <= simplify(_exp1))) then
        DEBUG(__FILE__, __LINE__, lprint(">> Case 6 @ averkov_lemma_7"));
        pos_coeff := convert(
            evalf(
                solve(
                    simplify(_exp3) = simplify(_exp1), alpha, 'maxsols'=1)
                 ),
            rational);
        N := ceil(1/2*subs(alpha=pos_coeff, simplify(_exp1)));
    end if;

    if not(assigned(pos_coeff)) and not(assigned(N)) then
        pos_coeff := 1;
        N := 1;
    end if;

$ifdef LOG_TIME
    END_LOG_TIME("averkov_lemma_7::compute_N_heuristic",5);
$endif

    if (N > N_guess) then
        g := AVERKOV_EXPR(N_guess);
        if checkPositivityOverSAS(semialgebraic_of_B, f - g, x) then
            N := N_guess;
            DEBUG(__FILE__, __LINE__, lprint(">> N_guess was ok"));
        else
            DEBUG(__FILE__, __LINE__, lprint(">> N_guess was not ok"));
        end if;
    end if;

$ifdef ENABLE_BINARY_SEARCH_AVKL
    g := AVERKOV_EXPR(N);
$else
$ifdef ENABLE_AVERKOV_CHECK
    g := AVERKOV_EXPR(N);
$endif
$endif
$ifdef ENABLE_AVERKOV_CHECK
    lprint(">> 1. Checking correctness of averkov_lemma_7", univariate_sign_set([B_poly,g-f],["nonnegative","nonnegative"],x));
$endif

$ifdef LOG_TIME
    START_LOG_TIME("averkov_lemma_7::compute_N_binary_search",6);
$endif
    #
    # Find the smallest N satisfying the lemma requirement
    # f - g > 0 over SemiAlgebraic(B)
    # we use a binary search to refine N
    #
$ifdef ENABLE_BINARY_SEARCH_AVKL
    DEBUG(__FILE__, __LINE__, lprint(">> N before ENABLE_BINARY_SEARCH_AVKL", evalf(N)));
    # Semantics
    #  - AVERKOV_EXPR(N_top) satisfies conditions
    #  - AVERKOV_EXPR(N_cur) will be checked to reduce exponent
    N_top := N;
    N_bot := 0;

    init_sos := isSOS(f - g, x);

    while N_top - N_bot > 1 do
        N_cur := iquo(N_top + N_bot, 2);
        DEBUG(__FILE__, __LINE__, lprint(">> Current N_top", N_top));
        DEBUG(__FILE__, __LINE__, lprint(">> Current N_bot", N_bot));
        DEBUG(__FILE__, __LINE__, lprint(">> Current N_cur", N_cur));

        g := AVERKOV_EXPR(N_cur);

        if init_sos then
            predicate_bin_search:= isSOS(f - g, x);
        else
            predicate_bin_search := checkPositivityOverSAS(semialgebraic_of_B, f - g, x);
        end if;

        if predicate_bin_search then
            N_top := N_cur;
        else
            N_bot := N_cur;
        end if;
    end do;

    if N_top = 0 and checkPositivityOverSAS(semialgebraic_of_B,f,x) then
        N := -1;
    else
        N := N_top;
    end if;

    DEBUG(__FILE__, __LINE__, lprint(">> N after ENABLE_BINARY_SEARCH_AVKL", evalf(N)));
$endif

    DEBUG(__FILE__, __LINE__, lprint(">> N: ", N));
$ifdef LOG_TIME
    END_LOG_TIME("averkov_lemma_7::compute_N_binary_search",6);
$endif

    if N = -1 then
$ifdef LOG_TIME
        END_LOG_TIME("averkov_lemma_7",0)
$endif
        return map(g_i -> 0, basis);
    else
$ifdef LOG_TIME
        END_LOG_TIME("averkov_lemma_7",0)
$endif
        return map(g_i -> AVERKOV_SOS(g_i, N), basis);
    end if;
end proc;
