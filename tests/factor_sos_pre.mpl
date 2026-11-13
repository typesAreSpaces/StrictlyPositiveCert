# Run from the repository root, with and without -DFACTOR_SOS_PRE,
# -DENABLE_STAGE_DIVISION and -DWEIFENG_OPTIMIZATION. No installed .mla
# is needed.
$ifndef EPS_FACTOR
$define EPS_FACTOR 17/20
$endif
$ifndef N_GUESS_AVKL
$define N_GUESS_AVKL 2
$endif
$ifndef N_GUESS_LS
$define N_GUESS_LS 2
$endif
$include "src/StrictlyPositiveCert.mpl";
# Access local helpers for direct regression coverage without exporting them.
kernelopts(opaquemodules=false):

checkCertificate := proc(f, basis, cert, x)
local h;
    if nops(cert) <> nops(basis) + 1 or
       expand(f - StrictlyPositiveCert:-dot_product([1, op(basis)], cert)) <> 0 then
        error "certificate does not reconstruct its input";
    end if;
    for h in cert do
        if SemiAlgebraic([h < 0], [x]) <> [] then
            error "certificate coefficient is not nonnegative: %1", h;
        end if;
    end do;
end proc:

# Controlled stages force different factors to be extracted at each boundary.
# Both bound-generator slots are used, so both coefficient lifts are checked.
Stages := module()
export spCertificates;
local dot_product := StrictlyPositiveCert:-dot_product;
local SplitNonnegativeFactor := StrictlyPositiveCert:-SplitNonnegativeFactor;
local A := x^2+1, B := x^2+2, C := x^2+3;
# Never shortcut the mocked stages; division is covered separately below.
local divisionCertificate := proc(poly, g, x)
    return false, 0, poly;
end proc;
local bound_poly := proc(basis, x)
    return [-x^2, 1, 1, 1, 2];
end proc;
local averkov_lemma_7 := proc(x, f, basis, g, N)
$ifdef FACTOR_SOS_PRE
    return [1, 0];
$else
    return [A, 0];
$endif
end proc;
local lower_bound_poly := proc(x, f, g)
$ifdef FACTOR_SOS_PRE
    if expand(f - (-x^4+C*(1-x^2))) <> 0 then
        error "first remainder was not reduced";
    end if;
    return x^2, 1/10;
$else
    return A*B*x^2, 1/10;
$endif
end proc;
local averkov_extended_lemma := proc(x, f, g, N, eps)
$ifdef WEIFENG_OPTIMIZATION
    if eps <> -1 then
        error "expected epsilon sentinel in Weifeng branch";
    end if;
$endif
$ifdef FACTOR_SOS_PRE
    if expand(f - (1-x^2)) <> 0 then
        error "final remainder was not reduced";
    end if;
    return 1;
$else
    return A*B*C;
$endif
end proc;
$include "src/spCertificates.mpl";
end module:

runTests := proc()
local A := x^2+1, B := x^2+2, C := x^2+3;
local basis, f, expected_multiplier, cert, test_case, S, value, H, eps;

for test_case in [[1, true], [-1, false], [0, false],
                  [(x^2+1)^3, true], [-(x^2+1)^3, false],
                  [(x^2+1)*(x^2-1), false], [x^2, false]] do
    if StrictlyPositiveCert:-isSOS(test_case[1], x) <> test_case[2] then
        error "incorrect strict positivity result for %1", test_case[1];
    end if;
end do;

for test_case in [
    [(x^2+1)*(1-x^2) + x^2+2, 1-x^2, true, x^2, 3],
    [3*(1-x^2), 1-x^2, true, 3, 0],
    [x^2*(1-x^4) + (x-1)^2, 1-x^4, true, x^2, (x-1)^2],
    [x^2+1, 1-x^4, true, 0, x^2+1],
    [x*(1-x^2) + 1, 1-x^2, false, 0, x*(1-x^2) + 1]
] do
    value, H, f := StrictlyPositiveCert:-divisionCertificate(
        test_case[1], test_case[2], x);
    if value <> test_case[3] or expand(H-test_case[4]) <> 0 or
       expand(f-test_case[5]) <> 0 then
        error "incorrect division certificate for %1", test_case[1];
    end if;
end do;
# quickCertificates must divide f by the generator, not the reverse.
basis := [1-x^2];
f := (x^2+1)*(1-x^2) + x^2+2;
value, cert := StrictlyPositiveCert:-quickCertificates(f, basis, x);
if not value then
    error "quick certificate not found";
end if;
checkCertificate(f, basis, cert, x);

S := SemiAlgebraic([x >= -2, x <= 2], [x]);
if not StrictlyPositiveCert:-checkPositivityOverSAS(S, (x^2+1)^3, x) or
   StrictlyPositiveCert:-checkPositivityOverSAS(S, (x^2+1)*(x^2-1), x) then
    error "incorrect positivity over an interval";
end if;
# All intervals still need their endpoints checked after roots are exhausted.
S := SemiAlgebraic([(x+2)*(x+1)*(x-1)*(x-2) <= 0], [x]);
if StrictlyPositiveCert:-checkPositivityOverSAS(S, -x*(x^2+1), x) then
    error "positivity check skipped a later interval";
end if;
S := SemiAlgebraic([x >= -2, x <= 2], [x]);
value := StrictlyPositiveCert:-computeMin(S, (x^2+1)^3+7, x);
if abs(value-8) > 1/10^6 then
    error "minimization changed the objective's values";
end if;
value := StrictlyPositiveCert:-computeMin(S, 7, x);
if value <> 7 then
    error "constant objective was not preserved";
end if;
value := StrictlyPositiveCert:-findEps(x,
    [1-x^2, (x^2+2)*(1-x^2)], (x^2+1)*(2-x^2));
if abs(evalf(value)-4) > 1/10^6 then
    error "epsilon computation changed generator values";
end if;
value := StrictlyPositiveCert:-find_constant_lower_bound_poly(
    (x^2+1)^2, 1, 1-x^2, x);
if abs(evalf(value)-999/1000) > 1/10^6 then
    error "lower-bound constant changed the original ratio";
end if;
# Each case: [f, g, f - H*g expected to be strictly positive over R].
# The first case used to divide by zero; the 1-1999/2000*x^2 case used to
# return a constant leaving a negative leading coefficient; the last one
# has f < 0 in the gap of S(g) and must fall back to positivity on S(g).
for test_case in [
    [x^3+2, 1-x^2, true],
    [x^3+3*x+5, 1-x^2, true],
    [(x^2+1)*(x+2), 1-x^2, true],
    [x^5-x+3, 1-x^2, true],
    [x^3+11/10, 1-x^2, true],
    [(x^2+1)^2*(x+2), 1-x^2, true],
    [1-1999/2000*x^2, 1-x^2, true],
    [(x^2-1)^3+1/1000-x^7/10^5, -(x^2-1)*(x^2-4), false]
] do
    H, eps := StrictlyPositiveCert:-lower_bound_poly(x, test_case[1], test_case[2]);
    f := expand(test_case[1] - test_case[2]*H);
    if SemiAlgebraic([H < 0], [x]) <> [] or
       not (type(degree(f, x), even) and lcoeff(f, x) > 0) or
       SemiAlgebraic([test_case[2] >= 0, f <= 0], [x]) <> [] then
        error "lower-bound construction is invalid for %1", test_case[1];
    end if;
    if StrictlyPositiveCert:-isSOS(f, x) <> test_case[3] then
        error "unexpected global positivity for %1", test_case[1];
    end if;
    if test_case[3] and eps <> -1 then
        error "eps_LS should be unused when f - H*g is SOS";
    end if;
end do;

basis := [-x^2/2, -x^2/2]:
$ifdef WEIFENG_OPTIMIZATION
f := A*(-x^2/2+B*C*(1-x^2)):
expected_multiplier := A*B*C:
$else
f := A*(-x^2/2+B*(-x^4+C*(1-x^2))):
expected_multiplier := A*B*x^2+A*B*C:
$endif
cert := Stages:-spCertificates(f, basis, x):
checkCertificate(f, basis, cert, x):
if expand(cert[1]-A*B*C) <> 0 or
   expand(cert[2]-(A+expected_multiplier)) <> 0 or
   expand(cert[3]-expected_multiplier) <> 0 then
    error "incorrect cumulative SOS factors";
end if:

# Exercise the real pipeline and its early returns as well.
for test_case in [
    [(x^2+1)*(x+2), [1-x^2]],
    [x+1, [-x^2]],
    [(x^2+1)*(x^2+2), [-x^2]],
    [(x^2+1)*(x+1), [-1]]
] do
    checkCertificate(test_case[1], test_case[2],
        StrictlyPositiveCert:-spCertificates(op(test_case), x), x);
end do:
$ifndef FACTOR_SOS_BENCHMARK
printf("SOS factor preprocessing tests passed\n"):
$endif
end proc:
$ifndef FACTOR_SOS_BENCHMARK
runTests():
quit:
$endif
