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

# Intersect fixed cells with all shifted generator nonnegative sets.
local epsilon_bad_cells := proc(cells0, gens, offset, x)
local cells := cells0, g;
    for g in gens do
        cells := univariate_intersect_cells(cells,
                    univariate_select_cells(univariate_sign_chart(g+offset,x),"nonnegative"));
        if cells=[] then break; end if;
    end do;
    return cells;
end proc;

# Return a positive dyadic eps such that f>0 whenever every g_i+2*eps>=0.
# With B supplied, only its compact nonnegative set needs to be considered.
# This is a safe separation value, not an optimal separation threshold.
local findEps := proc(x, gens, f, B := FAIL)
local bad_cells, eps := 1;
    bad_cells := univariate_select_cells(univariate_sign_chart(f,x),"nonpositive");
    if B<>FAIL then
        bad_cells := univariate_intersect_cells(bad_cells,
                        univariate_select_cells(univariate_sign_chart(B,x),"nonnegative"));
    end if;
    if epsilon_bad_cells(bad_cells,gens,0,x)<>[] then
        error "f must be strictly positive on the nonnegative set of the generators";
    end if;
    while epsilon_bad_cells(bad_cells,gens,2*eps,x)<>[] do
        eps := eps/2;
    end do;
    return eps;
end proc;
