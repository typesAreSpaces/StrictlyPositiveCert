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

local find_g_min := proc(basis, point, x)
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
local computeMinInterval := proc(poly, a, b, x)
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

local findEps := proc(basis, T, x)
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

# Obtains argmin_{g \in basis}(g(sample_point))
local gMinPoint := proc(x, basis, sample_point)
DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> basis", basis));
local g_min := basis[1];
local g;
  for g in basis do
    if (evalf(subs(x=sample_point, g-g_min)) <= 0) then
      g_min := g;
    end if;
  end do;
  return g_min;
end proc;

# Returns a point strictly contained in [left, right]
local samplePoint := proc(left, right)
  if (left = -infinity and right = infinity) then
    return 0;
  end if;
  if (left = -infinity) then
    return right - 1;
  end if;
  if (right = infinity) then
    return left + 1;
  end if;
  return (left + right)/2;
end proc;

# Outputs minimal eps such that
# f > 0 over Semialgebraic(basis + 2*eps)
local gMinSeq := proc(x, basis, f)
local i, j;
local l := nops(basis);
local sol;
local interval;
local _interval;

DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> basis @ gMinSeq", basis));
DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> f @ gMinSeq", f));

local partition_roots := {};
local num_roots := 0;
  for i from 1 to l-1 do
    for j from i+1 to l do
      for sol in map(_interval -> _interval[1], realroot(basis[i] - basis[j], 1/10000)) do
        if evalf(subs(x=sol, basis[i])<0) then
          partition_roots := partition_roots union {sol};
          num_roots := num_roots + 1;
        end if;
      end do;
    end do;
  end do;
  DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> partition_roots", evalf(partition_roots)));
  partition_roots := sort(convert(partition_roots, list));

  i := 1;
local min_epsilon := infinity, curr_epsilon, curr_g;
local S := SemiAlgebraic([-f>=0], [x]);
DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> f", f));
DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> S", S));
  for interval in map(_interval -> bound_info(x, _interval, 0), S) do
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> interval", evalf(interval)));
    while (i <= num_roots and evalf(partition_roots[i] <= interval[1])) do
      i := i + 1;
    end do;

local left_endpoint := interval[1];

    while (i <= num_roots and evalf(partition_roots[i] < interval[2])) do
      curr_g := gMinPoint(x, basis, samplePoint(left_endpoint, partition_roots[i]));
      curr_epsilon := -maximize(curr_g, x = left_endpoint .. partition_roots[i]);
      DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> curr_epsilon", evalf(curr_epsilon)));
      DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> curr_g", curr_g));
      if (evalf(min_epsilon > curr_epsilon)) then
        min_epsilon := curr_epsilon;
      end if;

      left_endpoint := partition_roots[i];
      i := i + 1;
    end do;

    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> sample_point", samplePoint(left_endpoint, interval[2])));
    curr_g := gMinPoint(x, basis, samplePoint(left_endpoint, interval[2]));
    curr_epsilon := -maximize(curr_g, x = left_endpoint .. interval[2]);
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> curr_epsilon", evalf(curr_epsilon)));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> curr_g", curr_g));
    if (evalf(min_epsilon > curr_epsilon)) then
      min_epsilon := curr_epsilon;
    end if;
  end do;

  DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> min_epsilon", evalf(min_epsilon)));

  return min_epsilon;
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
