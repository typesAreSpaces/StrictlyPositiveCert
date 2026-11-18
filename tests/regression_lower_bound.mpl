# Run from the repository root with the same -D values as the Makefile.
# Read source directly so this does not depend on a rebuilt .mla library.
read "src/StrictlyPositiveCert.mpl":
kernelopts(opaquemodules=false):

check_lower_bound := proc()
local x, scales, f, g, h, eps_LS, residual, certificates, multiplier;
    for scales in [[1/10^6, 1/1000], [1/10^12, 1/10^9], [1, 1]] do
        f := scales[1]*(x+2);
        g := scales[2]*(1-x^2);
        h, eps_LS := StrictlyPositiveCert:-lower_bound_poly(x, f, g);
        residual := expand(f-h*g);
        if not type(h, rational) or h <= 0 then
            error "expected a positive rational multiplier for scales %1", scales;
        end if;
        if SemiAlgebraic([g >= 0, residual <= 0], [x]) <> [] then
            error "lower-bound step lost positivity for scales %1", scales;
        end if;
        if not type(degree(residual, x), even) or lcoeff(residual, x) <= 0 then
            error "residual is not bounded below for scales %1", scales;
        end if;
    end do;

    # Retain the zero-G case on the singleton feasible set {0}.
    h, eps_LS := StrictlyPositiveCert:-lower_bound_poly(x, x+1, -x^2);
    if not type(h, rational) or h <= 0
       or SemiAlgebraic([-x^2 >= 0, x+1+h*x^2 <= 0], [x]) <> [] then
        error "zero-G case failed";
    end if;

    # Exercise the original failure through the public entry point.
    f := (x+2)/10^6;
    g := (1-x^2)/1000;
    certificates := StrictlyPositiveCert:-spCertificates(f, [g], x);
    if expand(f-StrictlyPositiveCert:-dot_product([1, g], certificates)) <> 0 then
        error "certificate identity failed";
    end if;
    for multiplier in certificates do
        if not type(multiplier, polynom(rational, x)) then
            error "certificate multiplier is not a rational polynomial";
        end if;
        if SemiAlgebraic([multiplier < 0], [x]) <> [] then
            error "certificate multiplier is not nonnegative";
        end if;
    end do;
    printf("PASS: three scaling cases, zero-G case, and public certificate regression\n");
end proc:
check_lower_bound():
