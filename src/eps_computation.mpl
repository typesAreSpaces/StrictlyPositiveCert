# Obtains argmin_{g \in basis}(g(sample_point))
local min_g_at_point := proc(x, basis, sample_point)
    DEBUG(__FILE__, __LINE__, lprint(">> basis", basis));
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
local findEps := proc(x, basis, f)
local i, j;
local l := nops(basis);
local sol;
local interval;
local _interval;

    DEBUG(__FILE__, __LINE__, lprint(">> basis @ findEps", basis));
    DEBUG(__FILE__, __LINE__, lprint(">> f @ findEps", f));

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
    DEBUG(__FILE__, __LINE__, lprint(">> partition_roots", evalf(partition_roots)));
    partition_roots := sort(convert(partition_roots, list));

    i := 1;
local min_epsilon := infinity, curr_epsilon, curr_g;
local S := SemiAlgebraic([-f>=0], [x]);
    DEBUG(__FILE__, __LINE__, lprint(">> f", f));
    DEBUG(__FILE__, __LINE__, lprint(">> S", S));
    for interval in map(_interval -> bound_info(x, _interval, 0), S) do
        DEBUG(__FILE__, __LINE__, lprint(">> interval", evalf(interval)));
        while (i <= num_roots and evalf(partition_roots[i] <= interval[1])) do
            i := i + 1;
        end do;

        local left_endpoint := interval[1];

        while (i <= num_roots and evalf(partition_roots[i] < interval[2])) do
            curr_g := min_g_at_point(x, basis, samplePoint(left_endpoint, partition_roots[i]));
            curr_epsilon := -maximize(curr_g, x = left_endpoint .. partition_roots[i]);
            DEBUG(__FILE__, __LINE__, lprint(">> curr_epsilon", evalf(curr_epsilon)));
            DEBUG(__FILE__, __LINE__, lprint(">> curr_g", curr_g));
            if (evalf(min_epsilon > curr_epsilon)) then
                min_epsilon := curr_epsilon;
            end if;

            left_endpoint := partition_roots[i];
            i := i + 1;
        end do;

        DEBUG(__FILE__, __LINE__, lprint(">> sample_point", samplePoint(left_endpoint, interval[2])));
        curr_g := min_g_at_point(x, basis, samplePoint(left_endpoint, interval[2]));
        curr_epsilon := -maximize(curr_g, x = left_endpoint .. interval[2]);
        DEBUG(__FILE__, __LINE__, lprint(">> curr_epsilon", evalf(curr_epsilon)));
        DEBUG(__FILE__, __LINE__, lprint(">> curr_g", curr_g));
        if (evalf(min_epsilon > curr_epsilon)) then
            min_epsilon := curr_epsilon;
        end if;
    end do;

    DEBUG(__FILE__, __LINE__, lprint(">> min_epsilon", evalf(min_epsilon)));

    return min_epsilon;
end proc;
