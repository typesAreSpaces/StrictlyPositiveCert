# We assume:
# 1. SemiAlgebraic(B_poly) is compact
# Return: list of sums of squares multipliers l
# such that f - dot_product(l, basis) > 0 over SemiAlgebraic(B_poly)
local averkov_lemma_7 := proc(x, f, basis, B_poly, N_guess)
$ifdef LOG_TIME
    INIT_START_LOG_TIME("averkov_lemma_7",0)
$endif
local _gamma, interval, lowerbound, upperbound;
local eps, tobe_disjoint_set;
local N, g, term;
local semialgebraic_of_B;
local R := PolynomialRing([x]);
local M, mu, m, N_list, temp_bound_N;
# DEBUG This is only to work out one particular example
local pos_coeff, _pos_coeff;
local _error := 1/1000;
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Start @averkov_lemma_7"));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> args"));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> f", f));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> basis", basis));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> B_poly", B_poly));

$ifdef LOG_TIME
    START_LOG_TIME("averkov_lemma_7::Minimization_f",1);
$endif
    semialgebraic_of_B := SemiAlgebraic(
        [B_poly >= 0], [x]);

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
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> M", M));
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
    # Find _gamma
    #
$ifdef LOG_TIME
    START_LOG_TIME("averkov_lemma_7::compute_gamma",2);
$endif
    _gamma := 1/2*max(
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
    _gamma := max(ceil(evalf(_gamma)), 1);

    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> _gamma", _gamma));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> B_poly", B_poly));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> f", f));

    #
    # Find exponent eps
    #
$ifdef LOG_TIME
    START_LOG_TIME("averkov_lemma_7::compute_eps",3);
$endif
local T := SemiAlgebraic([B_poly >= 0, f < 0], [x]);
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> basis", basis));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> f", f));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> semialgebraic_of_B", semialgebraic_of_B));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> T", evalf(T)));

    eps := 1/2*convert(evalf(gMinSeq(x, basis, f)), rational);

    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> eps", eps));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> eps", evalf(eps)));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> f", f));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> basis", basis));
$ifdef LOG_TIME
    END_LOG_TIME("averkov_lemma_7::compute_eps",3);
$endif

    #
    # Find mu
    #
$ifdef LOG_TIME
    START_LOG_TIME("averkov_lemma_7::compute_mu",4);
$endif
local semialgebraic_for_mu := SemiAlgebraic([B_poly >= 0, op(map(g_i -> g_i + EPS_FACTOR*eps >= 0, basis))], [x]);
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
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> mu", mu));
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
local _exp1 := (log(2*m*_gamma) - log(alpha*mu))/(log(_gamma + eps) - log(_gamma));
local _exp2 := (log(2*m*_gamma) - log(alpha*M))/(log(_gamma + eps) - log(_gamma));
local _exp3 := (log(alpha*M) - log(2*eps))/(log(_gamma + 2*eps) - log(_gamma + eps));
local pos_coeff1, pos_coeff2, N1, N2;
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> M", M));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> m", m));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> _gamma", _gamma));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> mu", mu));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> eps", evalf(eps)));
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
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> N pos_coeff", N, pos_coeff));
$ifdef LOG_TIME
    END_LOG_TIME("averkov_lemma_7::compute_N_heuristic",5);
$endif

    if (N > N_guess) then
        g := add(term,
                 term in map(g_i -> 1/pos_coeff*g_i*((g_i - _gamma)/(_gamma + eps))^(2*N_guess), basis));
        DEBUG(__FILE__, __LINE__,ENABLE_DEBUGGING, lprint(">> B_poly", B_poly));
        DEBUG(__FILE__, __LINE__,ENABLE_DEBUGGING, lprint(">> g - f", g - f));
        #if SemiAlgebraic([B_poly >= 0, g - f >= 0], [x]) = [] then
        if checkPositivityOverSAS(semialgebraic_of_B, f - g, x) then
            N := N_guess;
            DEBUG(__FILE__, __LINE__,ENABLE_DEBUGGING, lprint(">> N_guess was ok"));
        else
            DEBUG(__FILE__, __LINE__,ENABLE_DEBUGGING, lprint(">> N_guess was not ok"));
        end if;
    end if;

    g := add(term,
             term in map(g_i -> 1/pos_coeff*g_i*((g_i - _gamma)/(_gamma + eps))^(2*N), basis));
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
        local N_top := N;
        local N_bottom := 0;
        local N_old := N_top;
        local N_curr;

        while true do
            DEBUG(__FILE__, __LINE__,ENABLE_DEBUGGING, lprint(">> Current N_top", N_top));
            DEBUG(__FILE__, __LINE__,ENABLE_DEBUGGING, lprint(">> Current N_bottom", N_bottom));
            N_curr := iquo(N_top + N_bottom, 2);
            DEBUG(__FILE__, __LINE__,ENABLE_DEBUGGING, lprint(">> Current N_curr", N_curr));
            g := add(term,
                     term in map(g_i -> 1/pos_coeff*g_i*((g_i - _gamma)/(_gamma + eps))^(2*N_curr), basis));
            DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> B_poly", B_poly));
            DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> g", g));
            DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> f", f));
            DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> f - g", f-g));
            # TODO This is a bottle neck
            #if SemiAlgebraic([B_poly >= 0, g - f >= 0], [x]) = [] then
            if checkPositivityOverSAS(semialgebraic_of_B, f - g, x) then
                N_top := N_curr;
            else
                N_bottom := N_curr;
            end if;
            if N_curr = N_old then
                break;
            end if;
            N_old := N_curr;
        end do;

        if N_top = 0 and SemiAlgebraicSetTools:-IsEmpty([B_poly >= 0, f <= 0], R) then
            N := -1;
        else
            N := N_top;
        end if;

        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> N after ENABLE_BINARY_SEARCH_AVKL", evalf(N)));
    end if;

    # TODO Remove this, this is just for testing purposes
    N:=10;

    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> N: ", N));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> eps: ", eps));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> _gamma: ", _gamma));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> f: ", f));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> basis: ", basis));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> B_poly: ", B_poly));
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
        return map(g_i -> 1/pos_coeff*((g_i - _gamma)/(_gamma + eps))^(2*N), basis);
    end if;
end proc;

# We assume:
# 1. _poly is strictly positive over SemiAlgebraic([g >= 0], [x])
# 2. _poly has a lowerbound over \mathbb{R}, i.e., SemiAlgebraic([_poly <= 0], [x]) is bounded
# 3. SemiAlgebraic([g >= 0], [x]) is bounded
local averkov_extended_lemma := proc(x, _poly, g, N_guess, eps_LS)
$ifdef LOG_TIME
    INIT_START_LOG_TIME("averkov_extended_lemma",0)
$endif
local A;
local i, j;
local _gamma, eps;
local tobe_disjoint_set;
local N, N1, N2, poly := _poly, _g;
local pos_coeff := 1;
local semialgebraic_eps_lifted;
local m, mu, interval, lowerbound, upperbound;
local R := PolynomialRing([x]);

    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> poly", poly));
    # Check is poly is non-negative over \mathbb{R}
    #if SemiAlgebraic([poly < 0],[x]) = [] then
    if Isolate(poly) = [] then
        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Done because poly is a sos"));
$ifdef LOG_TIME
        END_LOG_TIME("averkov_extended_lemma",0)
$endif
        return 0;
    end if;

    # Since poly is not a non-negative
    # polynomial, we can assume the min value
    # for `poly` is negative

    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> poly", poly));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> g", g));

    if ENABLE_N_HEURISTIC then
        pos_coeff := findPositiveConstantAvoidExponent(poly, g);
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
    _gamma := convert(1/2*evalf(1.001*maximize(g)), rational);
    _gamma := max(_gamma, 1);
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> _gamma", _gamma));

    #
    # Compute exponent eps
    #
    if (eps_LS = -1) then
      DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> g:", g));
      DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> poly:", poly));
      eps := evalf(gMinSeq(x, [g], poly));
      DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> eps:", eps));
      eps := 1/2*convert(eps, rational);
    else
      eps := eps_LS;
    end;

    semialgebraic_eps_lifted := SemiAlgebraic(
        [g + EPS_FACTOR*eps >= 0], [x]);
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Done computation of semialgebraic_eps_lifted", semialgebraic_eps_lifted));

    #
    # Compute mu
    #
    mu := computeMin(semialgebraic_eps_lifted, poly, x);
    mu := convert(evalf(mu), rational);
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> mu", mu));

    #
    # Find exponent N
    #

    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Compute exponent N"));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> _gamma @ averkov_extended_lemma:", evalf(_gamma)));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> mu @ averkov_extended_lemma:", evalf(mu)));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> eps @ averkov_extended_lemma:", evalf(eps)));
#local _exp1 := (log(2*_gamma) - log(alpha*mu))/(log(_gamma + eps) - log(_gamma));
#local _exp2 := (log(-alpha*m) - log(2*eps))/(log(_gamma + 2*eps) - log(_gamma + eps));
    #pos_coeff := convert(
        #evalf(solve(_exp1 = _exp2, alpha, 'maxsols'=1)), rational);
    #N := ceil(1/2*subs(alpha=pos_coeff, _exp1));
local _exp1 := (log(2*_gamma) - log(mu))/(log(_gamma + eps) - log(_gamma));
    pos_coeff := 1;
    N := ceil(evalf(_exp1));

    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> N: ", N));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> N_guess: ", N_guess));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> pos_coeff: ", pos_coeff));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> _poly: ", _poly));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> g: ", g));

    if (N > N_guess) then
        _g := 1/pos_coeff*g*((g - _gamma)/(_gamma + eps))^(2*N_guess);
        # Check is _poly - _g is non-negative over \mathbb{R}
        #if SemiAlgebraic([_g - _poly >= 0], [x]) = [] then
        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Isolate(_poly - _g)", _poly - _g));
        if Isolate(_poly - _g) = [] then
          N := N_guess;
          DEBUG(__FILE__, __LINE__,ENABLE_DEBUGGING, lprint(">> N_guess was ok @ averkov_extended_lemma"));
        else
          DEBUG(__FILE__, __LINE__,ENABLE_DEBUGGING, lprint(">> N_guess was not ok @ averkov_extended_lemma"));
        end if;
    end if;

    while ENABLE_POST_OPT do
        _g := 1/pos_coeff*g*((g - _gamma)/(_gamma + eps))^(2*N);
        # Check is _poly - _g is non-negative over \mathbb{R}
        #if SemiAlgebraic([_g - _poly >= 0], [x]) = [] then
        if Isolate(_poly - _g) = [] then
          break;
        end if;
        N := N+1;
        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> N", N));
    end do;

    if ENABLE_BINARY_SEARCH_LS then
        local N_top := N;
        local N_bottom := 0;
        local N_old := N_top;
        local N_curr;

        while true do
            DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Current N_top", N_top));
            DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Current N_bottom", N_bottom));
            N_curr := iquo(N_top + N_bottom, 2);
            DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Current N_curr", N_curr));
            DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Current pos_coeff", pos_coeff));
            _g := 1/pos_coeff*g*((g - _gamma)/(_gamma + eps))^(2*N_curr);
            DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Current _g", _g));
            #if SemiAlgebraic([_g - _poly >= 0], [x]) = [] then
            if Isolate(_poly - _g) = [] then
                N_top := N_curr;
            else
                N_bottom := N_curr;
            end if;
            if N_curr = N_old then
                break;
            end if;
            N_old := N_curr;
        end do;

        if N_top = 0 and SemiAlgebraicSetTools:-IsEmpty([_poly <= 0], R) then
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
        return 1/pos_coeff*((g - _gamma)/(_gamma + eps))^(2*N);
    end if;
end proc;
