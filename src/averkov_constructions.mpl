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
local R := PolynomialRing([x]);
local M, mu, m;
# DEBUG This is only to work out one particular example
local pos_coeff;
local _error := 1/1000;
local T;
local semialgebraic_for_mu;
local _exp1, _exp2, _exp3;
local pos_coeff1, pos_coeff2, N1, N2;
local N_top, N_bot, N_cur;
local init_sos := false, predicate_bin_search;
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Start @averkov_lemma_7"));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> f", f));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> basis", basis));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> B_poly", B_poly));

$ifdef LOG_TIME
    START_LOG_TIME("averkov_lemma_7::Minimization_f",1);
$endif
    semialgebraic_of_B := SemiAlgebraic(
        [B_poly >= 0], [x]);
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> semialgebraic_of_B", semialgebraic_of_B));

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
    # DEBUG if problems
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> M", evalf(M)));
    # M := convert(evalf(M), rational);
    if evalf(M < 0) then
        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Done because f is strictly positive over SemiAlgebraic(B_poly)"));
$ifdef LOG_TIME
        END_LOG_TIME("averkov_lemma_7",0)
$endif
        return map(g_i -> 0, basis);
    end if;

    # The following is needed because in the `optimization step' (Find N) staring on
    # line 662 breaks otherwise
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
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> gamma", gamma));

    #
    # Find exponent eps
    #
$ifdef LOG_TIME
    START_LOG_TIME("averkov_lemma_7::compute_eps",3);
$endif
    T := SemiAlgebraic([B_poly >= 0, f < 0], [x]);
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> T", evalf(T)));

    eps := 1/2*convert(evalf(findEps(x, basis, f)), rational);
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> eps", evalf(eps)));
$ifdef LOG_TIME
    END_LOG_TIME("averkov_lemma_7::compute_eps",3);
$endif

    #
    # Find mu
    #
$ifdef LOG_TIME
    START_LOG_TIME("averkov_lemma_7::compute_mu",4);
$endif
    semialgebraic_for_mu := SemiAlgebraic([B_poly >= 0, op(map(g_i -> g_i + EPS_FACTOR*eps >= 0, basis))], [x]);
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> semialgebraic_for_mu", semialgebraic_for_mu));
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

    # DEBUG if problems
    mu := convert(evalf(mu), rational);
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> mu", mu));
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
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> _exp1", _exp1));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> _exp2", _exp2));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> _exp3", _exp3));

    # 1.
    if evalf(subs(alpha=1, simplify(_exp1) <= simplify(_exp2) and simplify(_exp2) <= simplify(_exp3))) then
        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Case 1 @ averkov_lemma_7"));
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
        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Case 2 @ averkov_lemma_7"));
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
        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Case 3 @ averkov_lemma_7"));
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
        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Case 4 @ averkov_lemma_7"));
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
        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Case 5 @ averkov_lemma_7"));
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
        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Case 6 @ averkov_lemma_7"));
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

    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> N before ENABLE_BINARY_SEARCH_AVKL", evalf(N)));
$ifdef LOG_TIME
    END_LOG_TIME("averkov_lemma_7::compute_N_heuristic",5);
$endif

    if (N > N_guess) then
        g := AVERKOV_EXPR(N_guess);
        #if SemiAlgebraic([B_poly >= 0, g - f >= 0], [x]) = [] then
        if checkPositivityOverSAS(semialgebraic_of_B, f - g, x) then
            N := N_guess;
            DEBUG(__FILE__, __LINE__,ENABLE_DEBUGGING, lprint(">> N_guess was ok"));
        else
            DEBUG(__FILE__, __LINE__,ENABLE_DEBUGGING, lprint(">> N_guess was not ok"));
        end if;
    end if;

    g := AVERKOV_EXPR(N);
    DEBUG(__FILE__, __LINE__, ENABLE_AVERKOV_CHECK, print(">> 1. Checking correctness of averkov_lemma_7", SemiAlgebraic([B_poly >= 0, g - f >= 0], [x])));

$ifdef LOG_TIME
    START_LOG_TIME("averkov_lemma_7::compute_N_binary_search",6);
$endif
    #
    # Find the smallest N satisfying the lemma requirement
    # f - g > 0 over SemiAlgebraic(B)
    # we use a binary search to refine N
    #
    if ENABLE_BINARY_SEARCH_AVKL then
        # Semantics
        #  - AVERKOV_EXPR(N_top) satisfies conditions
        #  - AVERKOV_EXPR(N_cur) will be checked to reduce exponent
        N_top := N;
        N_bot := 0;

        if isSOS(f - g) then
            init_sos := true;
        end if;

        while N_top - N_bot > 1 do
            N_cur := iquo(N_top + N_bot, 2);
            DEBUG(__FILE__, __LINE__,ENABLE_DEBUGGING, lprint(">> Current N_top", N_top));
            DEBUG(__FILE__, __LINE__,ENABLE_DEBUGGING, lprint(">> Current N_bot", N_bot));
            DEBUG(__FILE__, __LINE__,ENABLE_DEBUGGING, lprint(">> Current N_cur", N_cur));

            g := AVERKOV_EXPR(N_cur);

            if init_sos then
                predicate_bin_search:= isSOS(f - g);
            else
                #predicate_bin_search := SemiAlgebraic([B_poly >= 0, g - f >= 0], [x]) = [];
                predicate_bin_search := checkPositivityOverSAS(semialgebraic_of_B, f - g, x);
            end if;

            if predicate_bin_search then
                N_top := N_cur;
            else
                N_bot := N_cur;
            end if;
        end do;

        if N_top = 0 and SemiAlgebraicSetTools:-IsEmpty([B_poly >= 0, f <= 0], R) then
            N := -1;
        else
            N := N_top;
        end if;

        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> N after ENABLE_BINARY_SEARCH_AVKL", evalf(N)));
    end if;

    # FIX remove this
    #DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> hahaha"));
    #N := 10;

    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> N: ", N));
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

# We assume:
# 1. f is strictly positive over SemiAlgebraic([g >= 0], [x])
# 2. f has a lowerbound over \mathbb{R}, i.e., SemiAlgebraic([f <= 0], [x]) is bounded
# 3. SemiAlgebraic([g >= 0], [x]) is bounded
local averkov_extended_lemma := proc(x, f, g, N_guess, eps_LS)
$ifdef LOG_TIME
    INIT_START_LOG_TIME("averkov_extended_lemma",0)
$endif
local pos_coeff := 1;
local R := PolynomialRing([x]);
local averkov_poly;
local semialgebraic_eps_lifted;
local gamma, eps, mu, N;
local N_top, N_bot, N_cur;

    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> f", f));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> g", g));
    # Check is f is non-negative over \mathbb{R}
    #if SemiAlgebraic([f < 0],[x]) = [] then
    if isSOS(f) then
        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Done because f is a sos"));
$ifdef LOG_TIME
        END_LOG_TIME("averkov_extended_lemma",0)
$endif
        return 0;
    end if;

    # Since f is not a non-negative
    # polynomial, we can assume the min value
    # for `f` is negative

    if ENABLE_N_HEURISTIC then
        pos_coeff := findPositiveConstantAvoidExponent(f, g);
        if(pos_coeff > 0) then
            DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Done with pos_coeff"));
            DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> pos_coeff", pos_coeff));
$ifdef LOG_TIME
            END_LOG_TIME("averkov_extended_lemma",0)
$endif
            return pos_coeff;
        end if;
    end if;

    #
    # Compute gamma
    #
    # We just need a lowerbound, not the
    # tightest lowerbound [to discuss later]
    gamma := convert(1/2*evalf(1.001*maximize(g)), rational);
    gamma := max(gamma, 1);
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> gamma @ averkov_extended_lemma:", evalf(gamma)));

    #
    # Compute exponent eps
    #
    if (eps_LS = -1) then
        eps := evalf(findEps(x, [g], f));
        eps := 1/2*convert(eps, rational);
    else
        eps := eps_LS;
    end;
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> eps @ averkov_extended_lemma:", evalf(eps)));

    semialgebraic_eps_lifted := SemiAlgebraic(
        [g + EPS_FACTOR*eps >= 0], [x]);
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Done computation of semialgebraic_eps_lifted", semialgebraic_eps_lifted));

    #
    # Compute mu
    #
    mu := computeMin(semialgebraic_eps_lifted, f, x);
    mu := convert(evalf(mu), rational);
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> mu @ averkov_extended_lemma:", evalf(mu)));

    #
    # Find exponent N
    #

    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Compute exponent N"));
    N := ceil(evalf((log(2*gamma) - log(mu))/(log(gamma + eps) - log(gamma))));

    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> N: ", N));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> N_guess: ", N_guess));

    if (N > N_guess) then
        averkov_poly := AVERKOV_1_EXPR(g, N_guess);
        # Check is f - averkov_poly is non-negative over \mathbb{R}
        #if SemiAlgebraic([averkov_poly - f >= 0], [x]) = [] then
        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Isolate(f - averkov_poly)", f - averkov_poly));
        if isSOS(f - averkov_poly) then
            N := N_guess;
            DEBUG(__FILE__, __LINE__,ENABLE_DEBUGGING, lprint(">> N_guess was ok @ averkov_extended_lemma"));
        else
            DEBUG(__FILE__, __LINE__,ENABLE_DEBUGGING, lprint(">> N_guess was not ok @ averkov_extended_lemma"));
        end if;
    end if;

    while ENABLE_POST_OPT do
        averkov_poly := AVERKOV_1_EXPR(g, N);
        # Check is f - averkov_poly is non-negative over \mathbb{R}
        #if SemiAlgebraic([averkov_poly - f >= 0], [x]) = [] then
        if isSOS(f - averkov_poly) then
            break;
        end if;
        N := N+1;
        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> N", N));
    end do;

    if ENABLE_BINARY_SEARCH_LS then
        N_top := N;
        N_bot := 0;

        while N_top - N_bot > 1 do
            N_cur := iquo(N_top + N_bot, 2);
            DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Current N_top", N_top));
            DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Current N_bot", N_bot));
            DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Current N_cur", N_cur));
            averkov_poly := AVERKOV_1_EXPR(g, N_cur);
            #if SemiAlgebraic([averkov_poly - f >= 0], [x]) = [] then
            if isSOS(f - averkov_poly) then
                N_top := N_cur;
            else
                N_bot := N_cur;
            end if;
        end do;

        if N_top = 0 and SemiAlgebraicSetTools:-IsEmpty([f <= 0], R) then
            N := -1;
        else
            N := N_top;
        end if;

        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> N after ENABLE_BINARY_SEARCH", evalf(N)));
    end if;

    if N = -1 then
$ifdef LOG_TIME
        END_LOG_TIME("averkov_extended_lemma",0)
$endif
        return 0;
    else
$ifdef LOG_TIME
        END_LOG_TIME("averkov_extended_lemma",0)
$endif
        return AVERKOV_SOS(g, N);
    end if;
end proc;
