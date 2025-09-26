read "random_util.mpl";

strictlyLeftPolynomial := proc(x, isolated_points, _denominator)
local lowerbound := isolated_points[1];
local upperbound := isolated_points[nops(isolated_points)];
local rational_point := genRandomPoint(lowerbound, upperbound, _denominator);
    return (x - (lowerbound - rational_point));
end proc;

strictlyRightPolynomial := proc(x, isolated_points, _denominator)
local lowerbound := isolated_points[1];
local upperbound := isolated_points[nops(isolated_points)];
local rational_point := genRandomPoint(lowerbound, upperbound, _denominator);
    return -(x - (upperbound + rational_point));
end proc;

inbetweenPolynomial := proc(x, isolated_points, _denominator)
local num_isolated_points := nops(isolated_points);
local choice := 1 + (rand() mod (num_isolated_points - 1));
local a, b;
    a := genRandomPoint(isolated_points[choice], isolated_points[choice+1], _denominator);
    b := genRandomPoint(isolated_points[choice], isolated_points[choice+1], _denominator);
    while(evalb(a = b)) do
        a := genRandomPoint(isolated_points[choice], isolated_points[choice+1], _denominator);
        b := genRandomPoint(isolated_points[choice], isolated_points[choice+1], _denominator);
    end do;
    return (x-a)*(x-b);
end proc;

getArchimedeanPolynomial := proc(x, isolated_points, offset)
local max_point := max(abs(isolated_points[1]), abs(isolated_points[nops(isolated_points)])) + offset;
    return expand(max_point^2 - x^2);
end proc;
