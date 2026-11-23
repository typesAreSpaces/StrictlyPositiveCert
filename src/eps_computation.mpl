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

# Return a positive rational eps such that f>0 whenever every g_i+2*eps>=0.
# With B supplied, only its compact nonnegative set needs to be considered.
#
# Validity is monotone: a larger eps only enlarges {g_i>=-2*eps}. The valid
# values form an open interval (0,eps*), where eps* itself fails, so search for
# a large valid eps rather than the threshold: double while valid (or halve
# until valid), then bisect between the last valid and first invalid values.
# Larger eps makes the Averkov exponent converge faster. lo is always a value
# that passed the test, so the result is valid whenever the search stops.
local findEps := proc(x, gens, f, B := FAIL,
                      {bisection_steps::nonnegint := 3, max_doublings::nonnegint := 64})
local bad_cells, valid, lo, hi, mid, i;
    bad_cells := univariate_select_cells(univariate_sign_chart(f,x),"nonpositive");
    if B<>FAIL then
        bad_cells := univariate_intersect_cells(bad_cells,
                        univariate_select_cells(univariate_sign_chart(B,x),"nonnegative"));
    end if;
    # No point to separate: every eps is valid and doubling would not stop.
    if bad_cells=[] then return 1; end if;
    if epsilon_bad_cells(bad_cells,gens,0,x)<>[] then
        error "f must be strictly positive on the nonnegative set of the generators";
    end if;
    valid := e -> evalb(epsilon_bad_cells(bad_cells,gens,2*e,x)=[]);
    if valid(1) then
        # A nonempty bad set bounds eps*, so doubling terminates; the cap is a
        # safety net only.
        lo, hi := 1, 2;
        for i from 1 to max_doublings while valid(hi) do
            lo, hi := hi, 2*hi;
        end do;
        if valid(hi) then return hi; end if;
    else
        # Equivalent to bisecting (0,1): halve until the first valid value.
        hi := 1;
        while not valid(hi/2) do hi := hi/2; end do;
        lo := hi/2;
    end if;
    for i from 1 to bisection_steps do
        mid := (lo+hi)/2;
        if valid(mid) then lo := mid; else hi := mid; end if;
    end do;
    return lo;
end proc;
