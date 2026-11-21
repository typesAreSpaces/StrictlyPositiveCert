dot_product := proc(v1, v2)
$ifdef LOG_TIME
    INIT_START_LOG_TIME("dot_product",0)
$endif
local out := 0, i;
    for i from 1 to min(nops(v1), nops(v2)) do
        out := out + v1[i]*v2[i];
    end do;
$ifdef LOG_TIME
    END_LOG_TIME("dot_product",0)
$endif
    return out;
end proc;

# Check if poly is strictly positive
# over S
# S is a finite list of intervals
# poly is a polynomial
local checkPositivityOverSAS := proc(S, poly, x)
$ifdef LOG_TIME
    INIT_START_LOG_TIME("checkPositivityOverSAS",0)
$endif
local interval;
local local_poly := realroot(diff(poly, x), 1/10000);
    DEBUG(__FILE__, __LINE__, lprint(">> local_poly", evalf(local_poly)));
local curr_point;
local i, j := 1;
    for i from 1 to nops(S) do
        interval := bound_info(x, S[i], 0);
        DEBUG(__FILE__, __LINE__, lprint(">> Current interval", interval));
        if evalf(subs(x=interval[1], poly) <= 0) then
$ifdef LOG_TIME
            END_LOG_TIME("checkPositivityOverSAS",0)
$endif
            return false;
        end if;
        if evalf(subs(x=interval[2], poly) <= 0) then
$ifdef LOG_TIME
            END_LOG_TIME("checkPositivityOverSAS",0)
$endif
            return false;
        end if;

        if j > nops(local_poly) then
            break;
        end if;
        curr_point := (local_poly[j,1]+local_poly[j,2])/2;
        while j <= nops(local_poly) and evalf(interval[1] > curr_point) do
            curr_point := (local_poly[j,1]+local_poly[j,2])/2;
            j := j + 1;
        end do;
        DEBUG(__FILE__, __LINE__, lprint(">> curr_point", curr_point));
        while j <= nops(local_poly) and evalf(curr_point <= interval[2]) do
            curr_point := (local_poly[j,1]+local_poly[j,2])/2;
            DEBUG(__FILE__, __LINE__, lprint(">> Show evalf", evalf(subs(x=curr_point, poly)) ));
            if evalf(subs(x=curr_point, poly) <= 0) then
$ifdef LOG_TIME
                END_LOG_TIME("checkPositivityOverSAS",0)
$endif
                return false;
            end if;
            j := j + 1;
        end do;
    end do;
$ifdef LOG_TIME
    END_LOG_TIME("checkPositivityOverSAS",0)
$endif
    return true;
end proc;

# Despite its name, this checks strict positivity, not general SOS membership.
# A real polynomial without real roots has constant sign, so test x = 0 first.
isSOS := proc(poly, x)
    if not evalb(subs(x = 0, poly) > 0) then
        return false;
    end if;
    return evalb(Isolate(poly, x) = []);
end proc;

quickCertificates := proc(f, basis, x)
local certificates := map(gen -> 0, basis);
local q, r, i;
  for i from 1 to nops(basis) do
    if basis[i] = 0 then
      next;
    end if;
    # f = q*basis[i] + r gives the certificate multipliers directly.
    q := quo(f, basis[i], x);
    r := rem(f, basis[i], x);
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
