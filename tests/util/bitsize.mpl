#bitsize of integer
bitsizeI := proc(integer)
    if integer = 0 then
      return 1;
    end if;
    return ceil(evalf(log(abs(integer), 2))) + 1;
end proc;

# bitsize of rational number
bitsizeR := proc(rat)
local d := denom(rat);
    if (d = 1) then
        return bitsizeI(numer(rat));
    else
        #return max(bitsizeI(numer(rat)), bitsizeI(denom(rat)));
        return bitsizeI(numer(rat)) + bitsizeI(denom(rat));
    end if;
end proc;

# bitsize of polynomial with rational numbers
bitsizeP := proc(poly, x)
local _sum;
#return foldl((_x, _y) -> max(_x, _y), 0, op(map(bitsizeR, [coeffs(expanded)])));
    return add(_sum, _sum in map(bitsizeR, [coeffs(collect(poly, x))]));
end proc;
