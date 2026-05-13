with(SolveTools, SemiAlgebraic);
with(ListTools, FlattenOnce);

bound_info := proc(x, bound, eps)
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

Epsilon := proc(x, poly, basis, B_poly)
local _error := 1/1000;
local T := SemiAlgebraic([B_poly >= 0, poly < 0], [x]);
local rootsPositivePoly := [RealDomain:-solve(poly = 0, x)];
local lift_basis := lift -> map(_poly -> _poly - lift - _error, basis);
local eps_candidates :=
    FlattenOnce(map(
        proc(g_i)
            map(
                proc(bound)
                    interval := bound_info(x, bound, 0);
                    # TOCHECK
                    # This might introduce a bug if `lowerbound > upperbound`
                    # happens to be true for some reason
                    lowerbound := convert(evalf(interval[1]), rational);
                    upperbound := convert(evalf(interval[2]), rational);
                    simplify(maximize(g_i, x = lowerbound .. upperbound))
                end proc, T)
        end proc, basis));
local is_valid_eps :=
    proc(eps_candidate)
local i, j, check;
    # Check that for all roots
    # there is at least one negative point
    for i from 1 to nops(rootsPositivePoly) do
        # check is true if at least one poly evaluates to a negative number
        # after lifting with the eps_candidate and substitution at the ith root
        # of poly
        check := foldl(
            (x, y) -> x or y,
            false,
            op(map(poly -> subs(x = rootsPositivePoly[i], poly) < 0, lift_basis(eps_candidate))));
        if check = false then
            return false;
        end if;
    end do;
    return true;
end proc;
    return -1/2*min(map(
        proc(eps_candidate)
            if is_valid_eps(eps_candidate) then
                eps_candidate
            else
                infinity
            end if;
        end proc, select(_value -> _value < 0, eps_candidates)));
end proc;

# TestS
evalf(
    Epsilon(
        x,
        -81/16*x^4 + 369/16*x^3 - 63/2*x^2 + 19/2*x + 5,
        [x, x*(x - 1), (x - 1)*(x - 2), -x + 2],
        -2*x^2 + 4*x + 8));
g := -(x+2)*(x-4);
f1 := -(x+1)*(x-3);
nat := [x, x*(x - 1), (x - 1)*(x - 2), -(x - 2)];
nat2 := [x, x*(x - 1), (x - 1)*(x - 2), -10*(x - 2)];
nat3 := [10*x, x, x*(x - 1), (x - 1)*(x - 2), -(x - 2), -10*(x - 2)];

evalf(
    Epsilon(
        x,
        f1,
        nat,
        g));
evalf(
    Epsilon(
        x,
        f1,
        nat2,
        g));
evalf(
    Epsilon(
        x,
        f1,
        nat3,
        g));
