restart:
# Search the experimental archive before any installed copy of the package.
libname := "../StrictlyPositiveCert.mla", libname:
with(StrictlyPositiveCert);

interface(echo=0);
getenv(OMP_NUM_THREADS);
kernelopts(numcpus);

printlevel := -1;
#printlevel := 8;

# This file is intentionally separate from tests/run_regressions.py: the
# current implementation is expected to reproduce all four reported bugs.
# After fixes, each section should instead report PASS.
#
# An ERROR means the test itself could not establish its result (for example,
# an unexpected timeout). Such a result is not evidence that a bug was fixed.
# The script does not call restart or quit, and does not edit source/library
# files. It temporarily exposes private module members for finding 2 and
# restores the previous opaque-module setting in a finally clause.

run_review_findings := proc()
local x, f, basis, H, p, answer, value, good, redundant, item, witness;
local old_opaque := kernelopts(opaquemodules);
local bugs := 0, passes := 0, errors := 0;
local valid_certificate, report, report_error;

    # In one variable, a real polynomial is SOS iff it is nonnegative on R.
    # Check NONNEGATIVITY rather than calling StrictlyPositiveCert:-isSOS:
    # that helper actually checks strict positivity and rejects valid squares
    # with zeros. A reconstruction identity alone does not establish SOS.
    valid_certificate := proc(f, basis, H, x)
    local i, multiplier;
        if not type(H,list) or nops(H) <> nops(basis)+1 then return false; end if;
        for multiplier in H do
            if not type(multiplier,polynom(rational,x)) then return false; end if;
            if SolveTools:-SemiAlgebraic([multiplier<0],[x]) <> [] then
                return false;
            end if;
        end do;
        return evalb(expand(f-H[1]-add(H[i+1]*basis[i],i=1..nops(basis)))=0);
    end proc;

    report := proc(label, correct)
        if correct then
            passes := passes+1;
            printf("PASS: %s\n",label);
        else
            bugs := bugs+1;
            printf("BUG REPRODUCED: %s\n",label);
        end if;
    end proc;
    report_error := proc(label, details)
        errors := errors+1;
        printf("ERROR / INCONCLUSIVE: %s\n",label);
        lprint(details);
    end proc;

    kernelopts(opaquemodules=false);
    try
        # ------------------------------------------------------------------
        # 1. bound_poly pairs a positive EVEN-degree generator with an ODD one.
        #
        # The feasible set for [x^2,-x^5,x+1] is [-1,0], and f=x+2 >= 1 there.
        # Thus this is a valid strictly positive, compact/Archimedean example.
        #
        # The current code returns:
        #   [128*x^4-64*x^3+x+2, 64*x*(x-1)^2, 64, 0].
        # Its identity is correct, but the multiplier of x^2 is -256 at x=-1.
        # That multiplier is not SOS. The odd-pair construction should not
        # have used the redundant positive even generator x^2.
        #
        # Replacing -x^5 with -x^3 also exposes a negative power of x: the
        # current multiplier is 8*(x-1)^2/x, which is not even a polynomial.
        # Expected after a fix: BOTH calls return valid polynomial SOS lists.
        # ------------------------------------------------------------------
        printf("\n=== Finding 1: invalid SOS multipliers ===\n");
        try
            f := x+2;
            good := true;
            for basis in [[x^2,-x^5,x+1], [x^2,-x^3,x+1]] do
                H := timelimit(15,StrictlyPositiveCert:-spCertificates(
                    f,basis,x,max_search_time=5));
                lprint("basis",basis);
                lprint("returned certificate",H);
                lprint("identity residual (zero alone is insufficient)",
                       normal(f-StrictlyPositiveCert:-dot_product([1,op(basis)],H)));
                lprint("multiplier of x^2 at x=-1",subs(x=-1,H[2]));
                lprint("is that multiplier a rational polynomial?",
                       type(H[2],polynom(rational,x)));
                if not timelimit(10,valid_certificate(f,basis,H,x)) then
                    good := false;
                end if;
            end do;
            report("mixed-parity bound construction",good);
        catch:
            report_error("mixed-parity bound construction",[lastexception]);
        end try;

        # ------------------------------------------------------------------
        # 2. Interval positivity must not skip a minimum or a later interval.
        #
        # (a) p(1)>0 and p(3)>0, but p(2)=-1/100. The derivative has earlier
        # roots before [1,3]; the cursor logic skips the in-range minimum.
        # (b) p=-x has no derivative roots. It is positive on [-2,-1] but
        # negative on [1,2]. `break` skips the second interval altogether.
        # Expected after a fix: BOTH positivity decisions are false.
        # ------------------------------------------------------------------
        printf("\n=== Finding 2: false positivity decisions ===\n");
        try
            p := (x+1)^2*(x-2)^2-1/100;
            answer := timelimit(10,StrictlyPositiveCert:-univariate_positive_on_cells(
                p,[[1,3,true,true]],x));
            lprint("first positivity answer (expected false)",answer);
            lprint("exact negative witness p(2)",subs(x=2,p));
            good := evalb(answer=false);
            answer := timelimit(10,StrictlyPositiveCert:-univariate_positive_on_cells(
                -x,[[-2,-1,true,true],[1,2,true,true]],x));
            lprint("second positivity answer (expected false)",answer);
            lprint("exact negative witness at x=1",-1);
            if answer <> false then good := false; end if;
            report("interval and critical-point traversal",good);
        catch:
            report_error("interval and critical-point traversal",[lastexception]);
        end try;

        # ------------------------------------------------------------------
        # 3. Zero and positive constant generators crash bound preprocessing.
        #
        # Adding 0 or 1 to [1-x^2] does not change the feasible set [-1,1].
        # f=x+3/2 is strictly positive there. Both are valid redundant inputs.
        # Current errors: division by zero for 0; invalid coeff argument for 1.
        # Expected after a fix: all three calls return valid SOS certificates.
        # Preserve original basis indices when skipping redundant generators.
        # ------------------------------------------------------------------
        printf("\n=== Finding 3: redundant constant generators ===\n");
        try
            f := x+3/2;
            H := timelimit(15,StrictlyPositiveCert:-spCertificates(
                f,[1-x^2],x,max_search_time=5));
            if not timelimit(10,valid_certificate(f,[1-x^2],H,x)) then
                error "control certificate failed";
            end if;
            printf("Control without a redundant generator: valid\n");
            good := true;
            for redundant in [0,1] do
                basis := [redundant,1-x^2];
                try
                    H := timelimit(15,StrictlyPositiveCert:-spCertificates(
                        f,basis,x,max_search_time=5));
                    lprint("basis and certificate",basis,H);
                    if not timelimit(10,valid_certificate(f,basis,H,x)) then
                        good := false;
                    end if;
                catch:
                    item := [lastexception];
                    lprint("error for valid redundant generator",redundant,item);
                    # Recognize the specific preprocessing errors. An unrelated
                    # error/timeout is inconclusive, rather than this diagnosis.
                    value := sprintf("%a",item);
                    if StringTools:-Search("division by zero",value)=0 and
                    StringTools:-Search("coeff",value)=0 then
                        error "unexpected error in redundant-generator test: %1",item;
                    end if;
                    good := false;
                end try;
            end do;
            report("zero/nonnegative constant preprocessing",good);
        catch:
            report_error("zero/nonnegative constant preprocessing",[lastexception]);
        end try;

        # ------------------------------------------------------------------
        # 4. The direct monogenic API cannot handle f=x, g=-1.
        #
        # {g>=0} is empty, so strict positivity on that set is vacuous. The
        # input passes the API's even-degree/negative-leading-coefficient test.
        # An explicit certificate is:
        #   x = ((x+1)/2)^2 + ((x-1)/2)^2 * (-1).
        #
        # Current behavior: every product in the search family has degree 0;
        # the degree guard skips every k and reports exhausted search budgets.
        # Increasing max_k cannot fix this. A negative-constant shortcut is
        # needed in monogenicCertificates, as already exists in spCertificates.
        # Expected after a fix: the direct call returns a valid certificate.
        # ------------------------------------------------------------------
        printf("\n=== Finding 4: negative constant generator ===\n");
        try
            witness := [(x+1)^2/4,(x-1)^2/4];
            lprint("explicit witness",witness);
            lprint("witness is a valid certificate",valid_certificate(x,[-1],witness,x));
            good := true;
            try
                H := timelimit(15,StrictlyPositiveCert:-monogenicCertificates(
                    x,-1,x,max_search_time=5));
                lprint("direct monogenic result",H);
                good := timelimit(10,valid_certificate(x,[-1],H,x));
            catch:
                item := [lastexception];
                lprint("direct monogenic error",item);
                if StringTools:-Search("search exhausted",sprintf("%a",item))=0 then
                    error "unexpected error in negative-constant test: %1",item;
                end if;
                good := false;
            end try;
            report("negative-constant monogenic search",good);
        catch:
            report_error("negative-constant monogenic search",[lastexception]);
        end try;
    finally
        kernelopts(opaquemodules=old_opaque);
    end try;
    printf("\nSUMMARY: %d bug(s) reproduced; %d correct behavior(s); %d inconclusive error(s).\n",
           bugs,passes,errors);
    printf("Expected on the reviewed source: 4 bugs, 0 correct behaviors, 0 errors.\n");
end proc:

run_review_findings():
