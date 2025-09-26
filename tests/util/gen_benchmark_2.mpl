read "random_util.mpl";

basisNoBoundedPoly := proc(x, intervals, k)
local i, max_point, count_num_roots := 0;
local g1 := 1, g2 := -1;

    for i from 1 to nops(intervals) do
        # Isolated points
        if (nops(intervals[i]) = 1) then
            max_point := intervals[i][1];
            count_num_roots := count_num_roots + 1;
            g1 := g1*(x - intervals[i][1])^k;
            g2 := g2*(x - intervals[i][1])^k;
            # Regular intevals
        else
            max_point := intervals[i][2];
            count_num_roots := count_num_roots + 2;
            g1 := g1*(x - intervals[i][1])^k*(x - intervals[i][2])^k;
            g2 := g2*(x - intervals[i][1])^k*(x - intervals[i][2])^k;
        end if;
    end do;
    if evalb(count_num_roots mod 2 = 0) then
      g1 := g1*(x - (max_point + 1))^k;
      g2 := g2*(x - (max_point + 1))^k;
    end if;
    return [g1, g2];
end proc;

linearBasisNoBoundedPoly  := proc(x, intervals)
    return basisNoBoundedPoly(x, intervals, 1)
end proc;

