_pwd := currentdir();
currentdir(FileTools[ParentDirectory](__FILE__));

read "random_util.mpl";

uniformBasisGenerator := proc(x, intervals, k)
local i;
local g1 := 1, g2 := -1;

    for i from 1 to nops(intervals) do
        # Isolated points
        if (nops(intervals[i]) = 1) then
            #g2 := g2*(x - intervals[i][1])^(2*k);
            g2 := g2*(x - intervals[i][1])^(k+1);
            # Regular intevals
        else
            g1 := g1*(x - intervals[i][1])^k*(x - intervals[i][2])^k;
            g2 := g2*(x - intervals[i][1])^k*(x - intervals[i][2])^k;
        end if;
    end do;
    if evalb(g1 = 1) then
      return [x^2, g2];
    end if;
    return [g1, g2];
end proc;

linearBasisGenerator := proc(x, intervals)
    return uniformBasisGenerator(x, intervals, 1)
end proc;

currentdir(_pwd);
