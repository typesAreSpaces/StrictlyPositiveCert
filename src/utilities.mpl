dot_product := proc(v1, v2)
local out := 0, i;
$ifdef LOG_TIME
    INIT_START_LOG_TIME("dot_product")
    try
$endif
    for i from 1 to min(nops(v1), nops(v2)) do
        out := out + v1[i]*v2[i];
    end do;
    return out;
$ifdef LOG_TIME
    finally
        END_LOG_TIME("dot_product");
    end try;
$endif
end proc;

# Despite its name, this checks strict positivity, not general SOS membership.
# A real polynomial without real roots has constant sign, so test x = 0 first.
isSOS := proc(poly, x)
    if not evalb(subs(x = 0, poly) > 0) then
        return false;
    end if;
    return evalb(Isolate(poly, x) = []);
end proc;

quickCertificates := proc(f, gens, x)
local certificates := map(gen -> 0, gens);
local q, r, i;
    for i from 1 to nops(gens) do
        if gens[i] = 0 then
            next;
        end if;
        # f = q*gens[i] + r gives the certificate multipliers directly.
        q := quo(f, gens[i], x);
        r := rem(f, gens[i], x);
        if isSOS(q, x) and isSOS(r, x) then
            certificates[i] := q;
            return true, [r, op(certificates)];
        end if;
    end do;
    return false, [];
end proc;

# Split poly = fnonneg*fpos2 over the coefficient field. Factors without
# real roots have constant sign; normalize them to be strictly positive.
# Factors with real roots stay in fpos2, preserving rational coefficients.
SplitNonnegativeFactor := proc(poly, x)
local L, item, factor_poly, multiplicity, factor_sign;
local fnonneg := 1;
local fpos2;
    L := factors(poly);
    fpos2 := L[1];
    for item in L[2] do
        factor_poly, multiplicity := op(item);
        if sturm(factor_poly, x, -infinity, infinity) = 0 then
            factor_sign := sign(lcoeff(factor_poly, x));
            fnonneg := fnonneg*(factor_sign*factor_poly)^multiplicity;
            fpos2 := fpos2*factor_sign^multiplicity;
        else
            fpos2 := fpos2*factor_poly^multiplicity;
        end if;
    end do;
    return fnonneg, fpos2;
end proc;
