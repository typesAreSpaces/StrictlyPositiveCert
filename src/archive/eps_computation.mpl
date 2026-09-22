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

    return -7/10*convert(eps, rational);
end proc;
