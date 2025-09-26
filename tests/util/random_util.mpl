#restart;
#randomize();

with(RandomTools, Generate);

genRandomPoint := proc(lowerbound, upperbound, _denominator)
    return Generate(rational(range=lowerbound..upperbound,denominator=_denominator));
end proc;

intervalsGenerator := proc(x, num_isolated_points, lowerbound, upperbound, _denominator)
local i, _lowerbound := lowerbound, _uppperbound := upperbound;
local a, b;
    return [seq(
        if(evalb(modp(rand(), 2) = 0)) then
            a := genRandomPoint(_lowerbound, _uppperbound, _denominator);
            b := genRandomPoint(_lowerbound, _uppperbound, _denominator);
            _lowerbound := max(a, b);
            if(abs(_lowerbound - _uppperbound) >= 1/_denominator) then
                _uppperbound := _uppperbound + _lowerbound;
            end if;
            if(evalb(a = b)) then
                [a]
            else
                [min(a, b), max(a, b)]
            end if;
        else
            a := genRandomPoint(_lowerbound, _uppperbound, _denominator);
            _lowerbound := a;
            if(abs(_lowerbound - _uppperbound) >= 1/_denominator) then
                _uppperbound := _uppperbound + _lowerbound;
            end if;
            [a]
        end if,
        i=1..num_isolated_points)];
end proc;

getIsolatedPoints := proc(intervals)
local isolated_points := [];
local i;
    for i from 1 to nops(intervals) do
        isolated_points := [op(isolated_points), op(intervals[i])];
    end do;
    return isolated_points;
end proc;
