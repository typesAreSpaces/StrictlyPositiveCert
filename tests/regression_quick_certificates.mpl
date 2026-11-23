# Run from the repository root with -DENABLE_QUICK_CERT and the Makefile's -D values.
# Read source directly so this does not depend on a rebuilt .mla library.
read "src/StrictlyPositiveCert.mpl":

check_quick_certificates := proc()
local x, f, basis, ok, certificates, multiplier, test_case;
    for test_case in [
        [1-x^2, [2-x^4, 1/4-x^2], [3/4, 0, 1]],
        [3-x^4, [x, 0, 1-x^2], [2, 0, 0, x^2+1]]
                     ] do
        f, basis := test_case[1], test_case[2];
        ok, certificates := StrictlyPositiveCert:-quickCertificates(f, basis, x);
        if not ok or certificates <> test_case[3] then
            error "unexpected quick certificate for %1", test_case;
        end if;
        if nops(certificates) <> nops(basis)+1
        or expand(f-StrictlyPositiveCert:-dot_product([1, op(basis)], certificates)) <> 0 then
            error "quick certificate identity failed";
        end if;
        for multiplier in certificates do
            if not type(multiplier, polynom(rational, x))
            or SolveTools:-SemiAlgebraic([multiplier < 0], [x]) <> [] then
                error "invalid SOS multiplier";
            end if;
        end do;
        if StrictlyPositiveCert:-spCertificates(f, basis, x) <> certificates then
            error "public quick-certificate path returned unexpected output";
        end if;
    end do;

    for test_case in [[1-x^2, [x^2+1]], [-x^4-2, [1-x^2]], [x+1, [0]], [x+1, []]] do
        ok, certificates := StrictlyPositiveCert:-quickCertificates(test_case[1], test_case[2], x);
        if ok or certificates <> [] then
            error "expected quick-certificate rejection for %1", test_case;
        end if;
    end do;
    printf("PASS: quick certificate identities, SOS multipliers, public path, and rejection cases\n");
end proc:
check_quick_certificates():
