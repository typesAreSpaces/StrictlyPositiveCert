dot_product := proc(v1, v2)
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

# Check if poly is strictly positive
# over S
# S is a finite list of intervals
# poly is a polynomial
local checkPositivityOverSAS := proc(S, poly, x)
$ifdef LOG_TIME
    INIT_START_LOG_TIME("checkPositivityOverSAS",0)
$endif
local interval;
local local_poly := realroot(diff(poly, x), 1/10000);
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> local_poly", evalf(local_poly)));
local curr_point;
local i, j := 1;
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

# Return: If f is not a non-negative polynomial
# then return a positive constant c such that
# c f - g is non-negative
# Otherwise, return 0
local findPositiveConstantAvoidExponent := proc(f, g)
$ifdef LOG_TIME
    INIT_START_LOG_TIME("findPositiveConstantAvoidExponent",0)
$endif
local i, _args, curr_condition, conditions, pos_coeff;
local sol := solve(
    {c > 0, c * f - g >= 0},
    {x}, 'parametric', 'real', 'parameters' = {c});

    _args := op(sol);
    conditions := [];

    for i from 1 to nops(sol)/2 do
        if(evalb(_args[2*i] = [[x = x]])) then
            conditions := [evalf(_args[2*i - 1]), op(conditions)];
            pos_coeff := Minimize(c, map(`<=`@op, conditions))[1];
            if(pos_coeff = 0) then
$ifdef LOG_TIME
                END_LOG_TIME("findPositiveConstantAvoidExponent",0)
$endif
                return 1;
            else
                # We return the inverse because we actually need
                # to produce a multiplier for `g`
$ifdef LOG_TIME
                END_LOG_TIME("findPositiveConstantAvoidExponent",0)
$endif
                return 1/convert(pos_coeff, rational, exact);
            end if;
        end if;
        curr_condition := _args[2*i - 1];
        conditions := [op(0, curr_condition)(seq(map(v -> -evalf(v + 1/100), [op(curr_condition)]))),
                       op(conditions)];
    end do;

$ifdef LOG_TIME
    END_LOG_TIME("findPositiveConstantAvoidExponent",0)
$endif
    return 0;
end proc;

local isSOS := proc(poly)
    return evalb(Isolate(poly) = []);
end proc;
