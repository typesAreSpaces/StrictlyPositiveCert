$define ENABLE_DEBUGGING      true
$define ENABLE_VERIFICATION   false
$define ENABLE_BINARY_SEARCH  true
$define ENABLE_N_HEURISTIC    false
$define ENABLE_AVERKOV_CHECK  false
$define BINARY_SEARCH_TRIGGER 100
#$define LOG_TIME

$define DEBUG_EXIT lprint(">> Debugging, getting out"); return 0
$define DEBUG(F, L, y, x) if (y) then lprint(">> Debugging file ", F, " at line ", L); x; end if

$define START_LOG_TIME(X, S) stack_level:=stack_level+1;fd := FileTools:-Text:-Open("log_time.txt", append);local _log_time_S := time();FileTools:-Text:-WriteString(fd, cat("Start: ", X, " ", convert(stack_level, string), "\n"));FileTools:-Text:-Close(fd);
$define END_LOG_TIME(X, S) fd := FileTools:-Text:-Open("log_time.txt", append);FileTools:-Text:-WriteString(fd, cat("End: ", X, " ", convert(stack_level, string), "\nTime: ", convert(time() - _log_time_S, string), "\n"));FileTools:-Text:-Close(fd);stack_level:=stack_level-1;
$define INIT_START_LOG_TIME(X, S) local fd;START_LOG_TIME(X, S)

with(SolveTools, SemiAlgebraic);
with(RootFinding, Isolate);
with(Optimization, Maximize, Minimize);
with(ListTools, FlattenOnce);
with(RegularChains, SemiAlgebraicSetTools, PolynomialRing);

StrictlyPositiveCert := module() option package;

export spCertificates;

$ifdef LOG_TIME
local stack_level := -1;
$endif

# Returns a list of the form [g, h1, h2, index1, index2]
# where
# - g is the bound poly
# - h1, h2 are the sums of squares which make the linear combination for g
# - index1, index2 are the indexes in basis which make the linear combination for g
local bound_poly := proc(basis, x)
$ifdef LOG_TIME
    INIT_START_LOG_TIME("bound_poly",0)
$endif
local degrees, fst_coeffs, snd_coeffs, h1, h2;
local i, j;
    degrees := map(poly -> degree(poly, x), basis);
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> degress", degrees));

    fst_coeffs := map[indices](
        i -> coeff(basis[i], x^degrees[i]), basis);
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> fst_coeffs", fst_coeffs));

    snd_coeffs := map[indices](
        i ->
        if degrees[i] = 1 then
            subs(x = 0, basis[i])/abs(fst_coeffs[i])
        else
            coeff(basis[i]/abs(fst_coeffs[i]), x^(degrees[i] - 1))
        end if, basis);
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> snd_coeffs", snd_coeffs));

    for i from 1 to nops(basis) do
        if type(degrees[i], even) and fst_coeffs[i] < 0 then
$ifdef LOG_TIME
            END_LOG_TIME("bound_poly",0)
$endif
            return [basis[i], 1, 0, i, i];
        end if;
    end do;

    # At this point, every element in the basis has odd degree
    for i from 1 to nops(basis) - 1 do
        for j from i + 1 to nops(basis) do
            if fst_coeffs[i]*fst_coeffs[j] > 0 then
                next;
            end if;
            if degrees[i] = degrees[j] then
                h1 := x^2/abs(fst_coeffs[i]);
                h2 := (x + sign(fst_coeffs[i])*(1/2*snd_coeffs[i] + 1/2*snd_coeffs[j] + 1))^2
                /abs(fst_coeffs[j]);
            elif degrees[j] < degrees[i] then
                h1 := 1/abs(fst_coeffs[i]);
                h2 := x^(degrees[i] - degrees[j] - 2)
                *(x + sign(fst_coeffs[i])*(1/2*snd_coeffs[i] + 1/2*snd_coeffs[j] + 1))^2
                /abs(fst_coeffs[j]);
            else
                h1 := x^(degrees[j] - degrees[i] - 2)
                *(x - sign(fst_coeffs[i])*(1/2*snd_coeffs[i] + 1/2*snd_coeffs[j] + 1))^2
                /abs(fst_coeffs[i]);
                h2 := 1/abs(fst_coeffs[j]);
            end if;
$ifdef LOG_TIME
            END_LOG_TIME("bound_poly",0)
$endif
            return [expand(h1*basis[i] + h2*basis[j]), h1, h2, i, j];
        end do;
    end do;
end proc;

local dot_product := proc(v1, v2)
$ifdef LOG_TIME
    INIT_START_LOG_TIME("dot_product",0)
$endif
local out := 0, i;
    for i from 1 to min(nops(v1), nops(v2)) do
        out := out + v1[i]*v2[i];
    end do;
$ifdef LOG_TIME
    END_LOG_TIME("dot_product",0)
$endif
    return out;
end proc;

local bound_info := proc(x, bound, eps)
$ifdef LOG_TIME
    INIT_START_LOG_TIME("bound_info",0)
$endif
local i1, i2, j1, j2;
    # This is a bounded bounduality
    if(nops(bound) = 2) then
        i1 := simplify(op(bound[1])[1]);
        i2 := simplify(op(bound[1])[2]);
        j1 := simplify(op(bound[2])[1]);
        j2 := simplify(op(bound[2])[2]);
        if evalb(i1 = x) then
            if evalb(j1 = x) then
$ifdef LOG_TIME
                END_LOG_TIME("bound_info",0)
$endif
                return [min(i2, j2)+eps, max(i2, j2)-eps];
            else
$ifdef LOG_TIME
                END_LOG_TIME("bound_info",0)
$endif
                return [min(i2, j1)+eps, max(i2, j1)-eps];
            end if;
        else
            if evalb(j1 = x) then
$ifdef LOG_TIME
                END_LOG_TIME("bound_info",0)
$endif
                return [min(i1, j2)+eps, max(i1, j2)-eps];
            else
$ifdef LOG_TIME
                END_LOG_TIME("bound_info",0)
$endif
                return [min(i1, j1)+eps, max(i1, j1)-eps];
            end if;
        end if;
        # This is an equality or unbounded bounduality
    else
        i1 := simplify(op(bound[1])[1]);
        j1 := simplify(op(bound[1])[2]);
        if type(bound[1], `=`) then
            if evalb(i1 = x) then
$ifdef LOG_TIME
                END_LOG_TIME("bound_info",0)
$endif
                return [j1, j1];
            else
$ifdef LOG_TIME
                END_LOG_TIME("bound_info",0)
$endif
                return [i1, i1];
            end if;
        end if;
        if type(bound[1], `<=`) then
            if evalb(i1 = x) then
$ifdef LOG_TIME
                END_LOG_TIME("bound_info",0)
$endif
                return [-infinity, j1-eps];
            else
$ifdef LOG_TIME
                END_LOG_TIME("bound_info",0)
$endif
                return [i1+eps, infinity];
            end if;
        end if;
    end if;
end proc;

# Check if poly is strictly positive
# over S
# S is a finite list of intervals
# poly is a polynomial
local checkPositivityOverSAS := proc(S, poly, x)
$ifdef LOG_TIME
    INIT_START_LOG_TIME("checkPositivityOverSAS",0)
$endif
local local_poly := realroot(diff(poly, x), 1/10000);
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> local_poly", evalf(local_poly)));
local curr_point;
local i, j := 1;
local interval;
    for i from 1 to nops(S) do
        interval := bound_info(x, S[i], 0);
        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Current interval", interval));
        if evalf(subs(x=interval[1], poly) <= 0) then
$ifdef LOG_TIME
            END_LOG_TIME("checkPositivityOverSAS",0)
$endif
            return false;
        end if;
        if evalf(subs(x=interval[2], poly) <= 0) then
$ifdef LOG_TIME
            END_LOG_TIME("checkPositivityOverSAS",0)
$endif
            return false;
        end if;

        if j > nops(local_poly) then
            break;
        end if;
        curr_point := (local_poly[j,1]+local_poly[j,2])/2;
        while j <= nops(local_poly) and evalf(interval[1] > curr_point) do
            curr_point := (local_poly[j,1]+local_poly[j,2])/2;
            j := j + 1;
        end do;
        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> curr_point", curr_point));
        while j <= nops(local_poly) and evalf(curr_point <= interval[2]) do
            curr_point := (local_poly[j,1]+local_poly[j,2])/2;
            DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Show evalf", evalf(subs(x=curr_point, poly)) ));
            if evalf(subs(x=curr_point, poly) <= 0) then
$ifdef LOG_TIME
                END_LOG_TIME("checkPositivityOverSAS",0)
$endif
                return false;
            end if;
            j := j + 1;
        end do;
    end do;
$ifdef LOG_TIME
    END_LOG_TIME("checkPositivityOverSAS",0)
$endif
    return true;
end proc;

# Compute minimum of polynomial poly
# over semialgebraic set S
# S is a finite list of intervals
# poly is a polynomial
local computeMin := proc(S, poly, x)
$ifdef LOG_TIME
    INIT_START_LOG_TIME("computeMin",0)
$endif
local roots_poly := map(sol -> op(sol)[2], Isolate(diff(poly, x)));
local num_roots := nops(roots_poly);
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> poly", poly));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> roots_poly", roots_poly));
local curr_point;
local curr_min := infinity;
local i, j := 1;
local interval;
    for i from 1 to nops(S) do
        interval := bound_info(x, S[i], 0);
        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Current interval", evalf(interval)));

        curr_point := evalf(subs(x=convert(interval[1], rational), poly));
        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> curr_point", curr_point));
        if evalf(curr_point <= curr_min) then
            curr_min := curr_point;
        end if;
        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> curr_min", curr_min));

        curr_point := evalf(subs(x=convert(interval[2], rational), poly));
        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> curr_point", curr_point));
        if evalf(curr_point <= curr_min) then
            curr_min := curr_point;
        end if;
        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> curr_min", curr_min));

        while j <= num_roots and evalf(roots_poly[j] < interval[1]) do
            DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> j @1", j));
            j := j + 1;
        end do;

        while j <= num_roots and evalf(roots_poly[j] < interval[2]) do
            DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> j @2", j));
            curr_point := evalf(subs(x=convert(roots_poly[j], rational), poly));
            if evalf(curr_point <= curr_min) then
                curr_min := curr_point;
                DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> curr_arg_min", evalf(roots_poly[j])));
                DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> curr_min", evalf(curr_min)));
            end if;
            j := j + 1;
        end do;
    end do;
$ifdef LOG_TIME
    END_LOG_TIME("computeMin",0)
$endif
    return curr_min;
end proc;

local find_g_min;
local computeMinInterval;
export findEps;

    find_g_min := proc(basis, point, x)
    local i;
    local curr_index := 0, curr_min := infinity, curr;
        for i from 1 to nops(basis) do
            curr := subs(x = point, basis[i]);
            if evalf(curr < curr_min) then
                curr_min := curr;
                curr_index := i;
            end if;
        end do;
        return curr_index;
    end proc;

    # Compute the minimum value of polynomial `poly`
    # in closed interval [a, b]
    computeMinInterval := proc(poly, a, b, x)
    local roots_poly := map(sol -> op(sol)[2], Isolate(diff(poly, x)));
    local num_roots := nops(roots_poly);
    local curr_point, curr_min := infinity;
    local i := 1;

        curr_point := evalf(subs(x=a, poly));
        if evalf(curr_point <= curr_min) then
            curr_min := curr_point;
        end if;

        curr_point := evalf(subs(x=b, poly));
        if evalf(curr_point <= curr_min) then
            curr_min := curr_point;
        end if;

        while i <= num_roots and evalf(roots_poly[i] <= a) do
            i := i + 1;
        end do;

        while i <= num_roots and evalf(roots_poly[i] < b) do
            curr_point := evalf(subs(x=roots_poly[i], poly));
            if evalf(curr_point < curr_min) then
                curr_min := curr_point;
            end if;
            i := i + 1;
        end do;

        return curr_min;
    end proc;

    findEps := proc(basis, T, x)
    local i, j;
    local points := [], num_points;
    local interval;
    local left_end, right_end, g_min;
    local eps := -infinity, _eps;

        for i from 1 to nops(basis) - 1 do
            for j from i + 1 to nops(basis) do
                points :=
                [
                    op(points),
                    op(
                        select(_point -> evalf(subs(x = _point, basis[i]) < 0),
                               map(_isolated -> op(_isolated)[2],
                                   Isolate(basis[i] - basis[j], x)
                                  )
                              )
                      )
                ];
            end do;
        end do;

        points := ListTools:-MakeUnique(sort(points));
        num_points := nops(points);

        j := 1;
        for i from 1 to nops(T) do
            interval := bound_info(x, T[i], 0);

            left_end := interval[1];
            right_end := interval[2];

            while j <= num_points and evalf(points[j] <= left_end) do
                j := j + 1;
            end do;
            # At this point, j > num_points or left_end < points[j]

            while true do
                if j > num_points or evalf(points[j] >= right_end) then
                    g_min := basis[find_g_min(basis, (left_end + right_end)/2, x)];
                    _eps := -computeMinInterval(-g_min, left_end, right_end, x);
                    if eps < _eps then
                        eps := _eps;
                    end if;
                    break;
                else
                    g_min := basis[find_g_min(basis, (left_end + points[j])/2, x)];
                    _eps := -computeMinInterval(-g_min, left_end, points[j], x);
                    if eps < _eps then
                        eps := _eps;
                    end if;
                    left_end := points[j];
                    j := j + 1;
                end if;
            end do;
        end do;
        # TODO Figure out a `better` multiplier
        #
        return -7/10*convert(eps, rational);
    end proc;

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

    # |-
    # --------
    # old code
    # --------
    # M := -min(seq(minimize(poly, x = B[i][1] .. B[i][2]), i = 1 .. numelems(B)));
    # |-
    M := -computeMin(semialgebraic_of_B, f, x);
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
    _gamma := ceil(evalf(_gamma));

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
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> T", T));
    eps := 1/2*findEps(basis, T, x);
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> eps", eps));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> eps", evalf(eps)));
$ifdef LOG_TIME
    END_LOG_TIME("averkov_lemma_7::compute_eps",3);
$endif

$ifdef LOG_TIME
    START_LOG_TIME("averkov_lemma_7::compute_mu",4);
$endif
    # |-
    # --------
    # old code
    # --------
    # mu := min(seq(minimize(poly, x = temp[i][1] .. temp[i][2]), i = 1 .. numelems(temp)));
    # |-
local semialgebraic_for_mu := SemiAlgebraic([B_poly >= 0, op(map(g_i -> g_i + 17/10*eps >= 0, basis))], [x]);
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> semialgebraic_for_mu", semialgebraic_for_mu));
    mu := computeMin(semialgebraic_for_mu, f, x);

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

    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">>> N before ENABLE_BINARY_SEARCH", evalf(N)));
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
    # TODO Decide a threshold to trigger binary search optimization
    #if ENABLE_BINARY_SEARCH and N < BINARY_SEARCH_TRIGGER then
    if ENABLE_BINARY_SEARCH then
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

        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">>> N after ENABLE_BINARY_SEARCH", evalf(N)));
    end if;

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
# 1. SemiAlgebraic(g) is bounded
# Return: list of polynomials l
# such that poly - l[2] has a lowerbound
# over \mathbb{}
local Lower_bound_poly := proc(x, poly, g)
$ifdef LOG_TIME
    INIT_START_LOG_TIME("Lower_bound_poly",0)
$endif
local i;
local d_poly, c_poly;
local d_g, h;
local S, _point;
local G, C;
# The following is used in the
# minimization problem to find
# C in order to avoid the boundary
# points
local eps := 1/1000;

$ifdef LOG_TIME
    START_LOG_TIME("Lower_bound_poly::expand(poly)",1);
$endif
    d_poly := degree(expand(poly), x);
    c_poly := coeff(poly, x^d_poly);
$ifdef LOG_TIME
    END_LOG_TIME("Lower_bound_poly::expand(poly)",1);
$endif

    # If poly has a lowerbound over \mathbb{R}
    # then make no changes to poly
    if type(d_poly, even) and evalb(evalf(c_poly) > 0) then
        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> poly is bounded: "));
$ifdef LOG_TIME
        END_LOG_TIME("Lower_bound_poly",0)
$endif
        return 0;
    end if;

$ifdef LOG_TIME
    START_LOG_TIME("Lower_bound_poly::SemiAlgebraic(g)",2);
$endif
    S := map(
        bound -> bound_info(x, bound, eps),
        SemiAlgebraic([g >= 0], [x]));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> S", S));
    d_g := degree(expand(g), x);
    _point := S[1][1];
    # ToDiscuss Is it needed to make this point a rational number?
    _point := convert(evalf(_point), rational);
$ifdef LOG_TIME
    END_LOG_TIME("Lower_bound_poly::SemiAlgebraic(g)",2);
$endif

    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> d_g", d_g));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> g", g));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> d_poly", d_poly));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> poly", poly));

    # TODO Compute h using 'more diverse' _points
    if d_g <= d_poly then
        if type(d_poly - d_g, even) then
            h := (x - _point)^(d_poly - d_g + 2);
        else
            h := (x - _point)^(d_poly - d_g + 1);
        end if;
    else
        h := 1;
    end if;
    G := h*g;
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> h", h));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> G", G));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> To optimize", diff(poly,x)*G - poly*diff(G, x) ));

    #local opt_roots := [RealDomain:-solve(diff(poly,x)*G - poly*diff(G, x) = 0, x)];
    #local opt_roots := map(evalf, [RealDomain:-solve(diff(poly,x)*G - poly*diff(G, x) = 0, x)]);
#local opt_roots := map(out -> op(out)[2], evalf(Isolate(diff(poly,x)*G - poly*diff(G, x), maxprec=100, digits=30)));
    #Digits:=30;
local opt_roots := map(out -> op(out)[2], Isolate(diff(poly,x)*G - poly*diff(G, x), maxprec=1000, digits=30));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> opt_roots", opt_roots));

$ifdef LOG_TIME
    START_LOG_TIME("Lower_bound_poly::Minimization_problem",3);
$endif
    # We just need a lowerbound, not the
    # tightest lowerbound [to discuss later]
    C := 9/10*min(
        seq(
            if evalb(S[i][1] = S[i][2]) then
                # If we are minimizing over
                # an isolated point of S(g), any value of
                # C satisfy the minimization condition
                1
            else
                min(
                    map(
                        x_arg -> subs({x=x_arg}, poly/G),
                        select(_root -> evalf(S[i][1] <= _root) and evalf(_root <= S[i][2]), opt_roots)))
                #opt_roots))
                #minimize(
                #simplify(poly/G),
                #x = S[i][1] .. S[i][2])
            end if,
            i = 1 .. numelems(S)));
$ifdef LOG_TIME
    END_LOG_TIME("Lower_bound_poly::Minimization_problem",3);
$endif

    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> C", C));
    C := evalf(C);
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> C as float", C));
    C := convert(C, rational);
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> C as rational", C));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> h", h));
$ifdef LOG_TIME
    END_LOG_TIME("Lower_bound_poly",0)
$endif
    return C*h;
end proc;

# We assume:
# 1. _poly is strictly positive over SemiAlgebraic([g >= 0], [x])
# 2. _poly has a lowerbound over \mathbb{R}, i.e., SemiAlgebraic([_poly <= 0], [x]) is bounded
# 3. SemiAlgebraic([g >= 0], [x]) is bounded
local Last_step := proc(x, _poly, g)
$ifdef LOG_TIME
    INIT_START_LOG_TIME("Last_step",0)
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
    #if SemiAlgebraic([poly < 0],[x]) = [] then
    if Isolate(poly) = [] then
        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Done because poly is a sos"));
$ifdef LOG_TIME
        END_LOG_TIME("Last_step",0)
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
            END_LOG_TIME("Last_step",0)
$endif
            return pos_coeff;
        end if;
    end if;

    pos_coeff := 1;

    # We just need a lowerbound, not the
    # tightest lowerbound [to discuss later]
    _gamma := convert(1/2*evalf(1.001*maximize(g)), rational);

    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> _gamma", _gamma));

    #
    # Find exponent eps
    #
    # |-
    # --------
    # old code
    # --------
    #eps := -1/2*convert((Maximize(g, {poly <= 0})[1]), rational, exact);
    # |-
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Start computation of SemiAlg_poly", poly));
    #1 local SemiAlg_poly := SemiAlgebraic([poly<=0], [x]);
local _SemiAlg_poly := map(arg -> op(arg)[2], Isolate(poly));
local _old_point;
    # SemiAlg_poly is the semialgebraic set of -poly
    # Since poly has a lower bound in \mathbb{R}, then
    # SemiAlgebraic([poly <= 0], [x]) is bounded
local SemiAlg_poly := [];
    for i from 1 to nops(_SemiAlg_poly) do
        if type(i, even) then
            SemiAlg_poly := [op(SemiAlg_poly), [_old_point, _SemiAlg_poly[i]]];
        end if;
        _old_point := _SemiAlg_poly[i];
    end do;
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Done computation of SemiAlg_poly", SemiAlg_poly));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Start computation of eps"));
    eps := -1/2*max(
        map(proc(bound)
                #1 interval := bound_info(x, bound, 0);
                # TOCHECK
                # This might introduce a bug if `lowerbound > upperbound`
                # happens to be true for some reason
                #1lowerbound := convert(evalf(interval[1]), rational);
                #1upperbound := convert(evalf(interval[2]), rational);
                #1lowerbound := interval[1];
                #1upperbound := interval[2];
                lowerbound := bound[1];
                upperbound := bound[2];
                #print(lowerbound<=upperbound);
                #print(evalb(lowerbound<=upperbound));
                simplify(maximize(g, x = lowerbound .. upperbound))
            end proc,
            SemiAlg_poly));
    # DEBUG if problems
    eps := convert(evalf(eps), rational);
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> eps:", eps));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> g:", g));

    semialgebraic_eps_lifted := SemiAlgebraic(
        [g + 17/10*eps >= 0], [x]);
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Done computation of semialgebraic_eps_lifted", evalf(semialgebraic_eps_lifted)));
    #mu := min(
    #map(proc(bound)
    #interval := bound_info(x, bound, 0);
    ## TOCHECK
    ## This might introduce a bug if `lowerbound > upperbound`
    ## happens to be true for some reason
    #lowerbound := convert(evalf(interval[1]), rational);
    #upperbound := convert(evalf(interval[2]), rational);
    ##simplify(minimize(poly, x = lowerbound .. upperbound))
    #simplify(Minimize(poly, {lowerbound <= x, x <= upperbound})[1])
    #end proc,
    #semialgebraic_eps_lifted));
    mu := computeMin(semialgebraic_eps_lifted, poly, x);
    mu := convert(evalf(mu), rational);
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> mu", mu));
    # |-
    # --------
    # old code
    # --------
    #temp := Cap_set(S, [[-infinity, infinity]]);
    #mu := min(seq(minimize(poly, x = temp[i][1] .. temp[i][2]), i = 1 .. numelems(temp)));
    #m := ceil(evalf(minimize(poly))) - 1;
    #lprint(_gamma, eps, mu, m);
    #N := 1/2*max(
    #ceil((log(mu) - log(2*_gamma))/(log(_gamma) - log(_gamma + eps))),
    #ceil((log(-m) - log(2*eps))/(log(_gamma + 2*eps) - log(_gamma + eps))));
    # |-

    #m := ceil(evalf(minimize(poly))) - 1;
    m := ceil(min(map(point -> subs(x=op(point)[2], poly), Isolate(diff(poly, x))))) - 1;

    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> m", m));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> m", evalf(m)));

    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Compute exponent N"));
    #
    # Find exponent N
    #
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> _gamma", _gamma));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> mu", mu));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> eps", eps));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> m", m));
local _exp1 := (log(2*_gamma) - log(alpha*mu))/(log(_gamma + eps) - log(_gamma));
local _exp2 := (log(-alpha*m) - log(2*eps))/(log(_gamma + 2*eps) - log(_gamma + eps));
    pos_coeff := convert(
        evalf(solve(_exp1 = _exp2, alpha, 'maxsols'=1)), rational);
    N := ceil(1/2*subs(alpha=pos_coeff, _exp1));

    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> N: ", N));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> pos_coeff: ", pos_coeff));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> _poly: ", _poly));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> g: ", g));

    # TODO Decide a threshold to trigger binary search optimization
    if ENABLE_BINARY_SEARCH and N < BINARY_SEARCH_TRIGGER then
        #if ENABLE_BINARY_SEARCH then
        #if false then
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
            if SemiAlgebraic([_g - _poly >= 0], [x]) = [] then
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
        END_LOG_TIME("Last_step",0)
$endif
        return 0;
    else
$ifdef LOG_TIME
        END_LOG_TIME("Last_step",0)
$endif
        return 1/pos_coeff*((g - _gamma)/(_gamma + eps))^(2*N);
    end if;
end proc;

    spCertificates := proc(f, basis, x)
$ifdef LOG_TIME
        INIT_START_LOG_TIME("spCertificates",0)
$endif
        if SemiAlgebraic([f < 0], [x]) = [] then
            return [f, op(map(0, basis))];
        end if;
    local g, H2, f2, H3, f3, H4, certificates;

        g := bound_poly(basis, x);
        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Done with bound_poly"));

        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Poly f for averkov_lemma_7", f));
        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Poly g[1] for averkov_lemma_7", g[1]));
        H2 := averkov_lemma_7(x, f, basis, g[1], 70);
        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Done with averkov_lemma_7"));

        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> H2", H2));
        f2 := f - dot_product(basis, H2);
        DEBUG(__FILE__, __LINE__, ENABLE_AVERKOV_CHECK, print(">> 2. Checking correctness of averkov_lemma_7", SemiAlgebraic([g[1] >= 0, f2 <= 0], [x])));
        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> f2", f2));

        H3 := Lower_bound_poly(x, f2, g[1]);
        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Done with Lower_bound_poly"));

        f3 := f2 - g[1]*H3;
        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> f3", f3));
        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> g[1]", g[1]));
        H4 := Last_step(x, f3, g[1]);
        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> Done with Last_step"));

        certificates := H2;
        certificates[g[4]] := certificates[g[4]] + (H3+H4)*g[2];
        if g[3] <> 0 then
            certificates[g[5]] := certificates[g[5]] + (H3+H4)*g[3];
        end if;

        certificates := [f - dot_product(basis, certificates), op(certificates)];
        DEBUG(__FILE__, __LINE__, ENABLE_VERIFICATION, lprint(">> This should be zero", expand(f - dot_product([1, op(basis)], certificates))));
        DEBUG(__FILE__, __LINE__, ENABLE_VERIFICATION, lprint(">> Certificates found", op(certificates)));
$ifdef LOG_TIME
        END_LOG_TIME("spCertificates",0);
$endif
        return certificates;
    end proc;
end module;
