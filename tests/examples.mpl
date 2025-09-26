with(StrictlyPositiveCert, spCertificates);

interface(echo=0);
getenv(OMP_NUM_THREADS);
kernelopts(numcpus);

printlevel := -1;
#printlevel := 8;

read "util/check_weifeng.mpl";
read "util/check_realcertify.mpl";

# Load implementation

# Load RealCertify
currentdir(homedir);
currentdir("Documents/GithubProjects/RealCertify");
read "multivsos/multivsos.mm";

# Original Weifeng examples
weifengExamples := proc()
    weifengExamples1();
    weifengExamples2();
    weifengExamples3();
end proc;

# Original Weifeng examples
weifengExamples1 := proc()
local x, f, g_1, g_2, G;

    f := -x*(x - 1)^3*(x - 2)^2;
    g_1 := x*(x - 1/2)*(x - 1)^2*(x - 2);
    g_2 := -x*(x - 1)*(x - 2);
    G := [g_1, g_2];

    try
        checkRealCertify(f, G, "Test 1");
    catch:
        printf(">> RealCertify fails Test 1\n");
    end try;
    try
        checkWeifeng(f, G, x, "Test 1");
    catch:
        printf(">> Weifeng approach fails Test 1\n");
    end try;
end proc;

weifengExamples2 := proc()
local x, f, g_1, g_2, G;
    f := -x*(x - 3);
    g_1 := x*(x - 1)*(x - 2)*(x - 3);
    g_2 := -x*(x - 1)*(x - 2)*(x - 3);
    G := [g_1, g_2];

    try
        checkRealCertify(f, G, "Test 2");
    catch:
        printf(">> RealCertify fails Test 2\n");
    end try;
    try
        checkWeifeng(f, G, x, "Test 2");
    catch:
        printf(">> Weifeng approach fails Test 2\n");
    end try;
end proc;

weifengExamples3 := proc()
local x, f, g_1, g_2, G;
    f := -x + 10;
    g_1 := (x - 2)^3;
    g_2 := -(x - 2)^3;
    G := [g_1, g_2];
    try
        checkRealCertify(f, G, "Test 3");
    catch:
        printf(">> RealCertify fails Test 3\n");
    end try;
    try
        checkWeifeng(f, G, x, "Test 3");
    catch:
        printf(">> Weifeng approach fails Test 3\n");
    end try;
end proc;

# RealCertify struggles with this example
# In general, we conjecture that RealCertify requires
# an Archimedean polynomial is the basis given consists only of
# polynomial whose semialgebraic set is not bounded
# But our approach succeeds
realCertifyIssues1 := proc()
local x;
    try
        checkRealCertify(x+1, [x*(x-1)*(x-2), -x*(x-1)*(x-2)], "Test");
    catch:
        printf(">> RealCertify fails with test\n");
    end try;

    checkWeifeng(x+1, [x*(x-1)*(x-2), -x*(x-1)*(x-2)], x, "Test");
end proc;

oldProblematicExamples1 := proc()
    checkRealCertify(
        x+101/100,
        [100-x^2,(x-118/25)*(x-657/100)*(x-89/10)*(x-413/20)*(x-561/25),
         (-x+118/25)*(x-657/100)*(x-89/10)*(x-413/20)*(x-561/25)], "Works");
    checkWeifeng(
        x+101/100,
        [(x-118/25)*(x-657/100)*(x-89/10)*(x-413/20)*(x-561/25),
         (-x+118/25)*(x-657/100)*(x-89/10)*(x-413/20)*(x-561/25)], x, "Works");
    checkRealCertify(
        x+53/100,
        [100-x^2,(x-577/100)*(x-317/50)*(x-1069/100)*(x-1443/100)*(x-381/25),
         -(x-577/100)*(x-317/50)*(x-1069/100)*(x-1443/100)*(x-381/25)], "Works");
    checkWeifeng(
        x+53/100,
        [(x-577/100)*(x-317/50)*(x-1069/100)*(x-1443/100)*(x-381/25),
         -(x-577/100)*(x-317/50)*(x-1069/100)*(x-1443/100)*(x-381/25)], x, "Now works!");
end proc;

oldProblematicExamples2 := proc()
    checkWeifeng(
        (x-7)*(x - (7 + 1/10)),
        [(x-118/25)*(x-657/100)*(x-89/10)*(x-413/20)*(x-561/25),
         (-x+118/25)*(x-657/100)*(x-89/10)*(x-413/20)*(x-561/25)],
        x, "PASSED");

    checkWeifeng(
        (x-7)*(x - (7 + 1/10)),
        [(x-118/25)*(x-657/100)*(x-89/10)*(x-413/20)*(x-561/25),
         (-x+118/25)^3*(x-657/100)*(x-89/10)*(x-413/20)*(x-561/25)],
        x, "FAILS");

    checkWeifeng(
        (x-10)*(x - 15),
        [(x-118/25)*(x-657/100)*(x-89/10)*(x-413/20)*(x-561/25),
         (-x+118/25)^3*(x-657/100)*(x-89/10)*(x-413/20)*(x-561/25)],
        x, "PASSED");
end proc;

realcertify_problematic_example := proc()
    checkRealCertify(
        -x^2+845306/625,
        [-x^2+839056/625,(x-61/25)*(x-483/50)*(x-1013/100)*(x-1539/100)*(x-1027/50)*(x-866/25)*(x-891/25), -(x-61/25)*(x-483/50)*(x-1013/100)*(x-1539/100)*(x-1027/50)*(x-866/25)*(x-891/25)],
#[-x^2+839056/625],
        x, "??");
end proc;

natGenExamples := proc()
    checkWeifeng((x+1)*x*(x-1), [x, x*(x-1), (x-1)*(x-2), -(x-2)], x, "A natural generator");
end proc;

complex_examples := proc()
    # -----------------------------------------------
    # Description:
    # Strictly positive polynomial over
    # saturated quadratic module using
    # natural generator basis

    #checkWeifeng(
    #-81/16*x^4 + 369/16*x^3 - 63/2*x^2 + 19/2*x + 5,
    #[x, x*(x-1), (x-1)*(x-2), -(x-2), 10 - x^2],
    #x, "PASSED");

    #checkWeifeng(
    #-81/16*x^4 + 369/16*x^3 - 63/2*x^2 + 19/2*x + 5,
    #[x, x*(x-1), (x-1)*(x-2), -(x-2)],
    #x, "PASSED");

    # RealCertify cannot handle this example
    #checkRealCertify(
    #-81/16*x^4 + 369/16*x^3 - 63/2*x^2 + 19/2*x + 5,
    #[x, x*(x-1), (x-1)*(x-2), -(x-2), 10 - x^2],
    #x, "FAILED");
    # -----------------------------------------------

    # -----------------------------------------------
    #checkWeifeng(
    #-(x+3)*(x+2)*(x+1)*(x-1)*(x-2)*(x-3),
    #[(x+3), (x+3)*x^2*(x-3), -(x-3)],
    #x, "FAILED");
    # -----------------------------------------------

    # -------------------------------------------------------------
    # Description:
    # Computing certificates of Archimedean polynomial
    # over non-saturated basis using generalized natural generators

    #checkWeifeng(
    #100-x^2,
    #[(x+3), (x+3)*x^2*(x-3), -(x-3)],
    #x, "PASSED");
    #checkWeifeng(
    #55-x^2,
    #[(x+3), (x+3)*x^2*(x-3), -(x-3)],
    #x, "PASSED");
    #checkWeifeng(
    #53-x^2,
    #[(x+3), (x+3)*x^2*(x-3), -(x-3)],
    #x, "TODO");
    # -------------------------------------------------------------

    # -------------------------------------------------------------
    # Description:
    # Computing certificates of strictly positive polynomial
    # over non-saturated basis using generalized natural generators

    #checkWeifeng(
    #(x+2)*(x+1),
    #[(x+3), (x+3)*x^2*(x-3), -(x-3)],
    #x, "PASSED");
    #checkRealCertify(
    #(x+2)*(x+1),
    #[(x+3), (x+3)*x^2*(x-3), -(x-3)],
    #x, "FAILED");
    #checkRealCertify(
    #(x+2)*(x+1),
    #[(x+3), (x+3)*x^2*(x-3), -(x-3), 10 - x^2],
    #x, "PASSED");
    # -------------------------------------------------------------

    # -------------------------------------------------------------
    # Description:
    # Computing certificates of products of generators
    # over non-saturated basis using generalized natural generators

    #checkWeifeng(
    #-(x+3)^3*(x-3)^3,
    #[(x+3)^3, (x+3)*x^2*(x-3), -(x-3)^3],
    #x, "PASSED");

    #checkWeifeng(
    #-(x+3)^2*x^2*(x-3)^2,
    #[(x+3), (x+3)*x^2*(x-3), -(x-3)],
    #x, "TODO");

    #checkWeifeng(
    #-(x+3)^2*(x-3)^2,
    #[(x+3), (x+3)*(x-3), -(x-3)],
    #x, "PASSED");

    #checkWeifeng(
    #-(x+3)^6*(x-3)^6,
    #map(poly -> 1/1*poly, [(x+3)^5, (x+3)^5*(x-3)^5, -(x-3)^5]),
    #x, "TODO");

    #checkWeifeng(
    #14000 - x^2,
    #map(poly -> 1/1*poly, [(x+3)^5, (x+3)^5*(x-3)^5, -(x-3)^5]),
    #x, "PASSED");

    #checkRealCertify(
    #-(x+(3+962/100000))*(x-(3+962/100000)),
    #[(x+3)^3, (x+3)^3*(x-3)^3, -(x-3)^3],
    #x, "FAILED");

    #map(k ->
    #checkWeifeng(
    #1000 - x^2,
    #[(x+3)^k, (x+3)^3*(x-3)^k, -(x-3)^k],
    #x, "TODO")
    #, [1, 3]);

    #checkRealCertify(
    #-(x+(3+962/100000))*(x-(3+962/100000)),
    #[(x+3)^3, (x+3)^3*(x-3)^3, -(x-3)^3, 9 - x^2],
    #x, "PASSED");

    #checkWeifeng(
    #-(x+3)*(x-3),
    #[(x+3), (x+3)*x^2*(x-3), -(x-3)],
    #x, "PASSED");
    #checkWeifeng(
    #(x+2)*(x+1)*(x-1)*(x-2),
    #[(x+3), (x+3)*x^2*(x-3), -(x-3)],
    #x, "PASSED");
    #checkRealCertify(
    #(x+2)*(x+1)*(x-1)*(x-2),
    #[(x+3), (x+3)*x^2*(x-3), -(x-3)],
    #x, "PASSED");
    # -------------------------------------------------------------
end proc;

close_to_semialgebraic := proc(eps)
    # -------------------------------------------------------------
    # Description:
    # Computing certificates of strictly positive polynomials
    # that are close to the semialgebraic set of the basis
local i;

printf("\n>> Start RealCertify\n");
    for i from 0 to 4 do
        try
            checkRealCertify(
                1+eps+x,
                [(1-x^2)^(2*i+1), 4 - x^2],
                "PASSED");
        catch:
            printf(">> RealCertify fails test with k: %d\n", 2*i+1);
        end try;
    end do;

printf("\n>> Start Weifeng\n");
    for i from 0 to 4 do
        try
            checkWeifeng(
                1+eps+x,
                [(1-x^2)^(2*i+1)],
                x, "PASSED");
        catch:
            printf(">> Weifeng fails test with k: %d\n", 2*i+1);
        end try;
    end do;
end proc;

simple_test1 := proc()
    checkWeifeng(
        (x+3)*(x+2),
        [-(x+5)*(x+4)*(x+3)*(x+2)*(x+1)*(x-1)*(x-2)*(x-3)*(x-4)*(x-5), (x+5)*(x+4)*(x+3)*(x+2)*(x+1)*(x-1)*(x-2)*(x-3)*(x-4)*(x-5)],
        x, "PASSED");
end proc;

simple_test2 := proc()
    checkWeifeng(
        -(x-1)^2,
        [x-1, -(x-1)],
        x, "PASSED");

    checkWeifeng(
        x+2,
        [x-1, -(x-1)],
        x, "PASSED");

    checkWeifeng(
        x+100,
        [x-1, -(x-1)],
        x, "PASSED");
end proc;

example_11_paper := proc()
    checkWeifeng(
        -26*x^7 + 13*x^6 + 87*x^5 + 49*x^4 - 464*x^3 + 1512*x^2 - 2211*x + 1092,
        [x*(x-1/2)*(x-1)^2*(x-2), -x*(x-1)*(x-2)],
        x, "PASSED");
end proc;

problematic_case := proc() 
    checkRealCertify(
        -x^2+4848041/10000,
        [(x-241/25)*(x-329/20)*(x-478/25)*(x-1979/100)*(x-2079/100), (-x+241/25)*(x-329/20)*(x-478/25)*(x-1979/100)*(x-2079/100)],
        x, "??");
    checkWeifeng(
        -x^2+4848041/10000,
        [(x-241/25)*(x-329/20)*(x-478/25)*(x-1979/100)*(x-2079/100), (-x+241/25)*(x-329/20)*(x-478/25)*(x-1979/100)*(x-2079/100)],
        x, "??");
end proc;

comparison_realcertify_fails := proc()
local f:=x+171/10;
#local g:=[(x-44/5)*(x-1289/100)*(x-2403/100)*(x-3099/100)*(x-3199/100), (-x+44/5)*(x-1289/100)*(x-2403/100)*(x-3099/100)*(x-3199/100), 10893401/10000-x^2];
local g:=[(x-44/5)*(x-1289/100)*(x-2403/100)*(x-3099/100)*(x-3199/100), (-x+44/5)*(x-1289/100)*(x-2403/100)*(x-3099/100)*(x-3199/100)];
    checkRealCertify(f, g,  x, "??");
    checkRealCertify(f, [op(g), (3199/100)^2 + 1 - x^2],  x, "??");
    checkWeifeng(f, g, x, "??");
end proc;

# ----------
# Examples
# ----------

printf("\n>> Start Examples\n");

#weifengExamples();
#weifengExamples1();
#weifengExamples2();
#weifengExamples3();

#realCertifyIssues1();
#simple_test1();
#simple_test2();
#example_11_paper();
#problematic_case();
#realcertify_problematic_example();
#natGenExamples();
#complex_examples();
#close_to_semialgebraic(1/2);
#close_to_semialgebraic(1/3);
#close_to_semialgebraic(1/5);
#comparison_realcertify_fails();

#checkWeifeng(x+2, [-(x+1), (x-1)], x, "PASSED");
#checkWeifeng(x+1/2, [-(x+1), (x-1)], x, "PASSED");

#checkWeifeng(-(x-5), [x-1, -(x+2)*(x+1)*(x-2)*(x-3), (x+1)*(x-4)], x, "PASSED");

#g1 := x + 3;
#g2 := (x + 2)*(x + 1);
#g3 := (x - 2)*(x - 1);
#g4 := -(x - 3);
#sigma1 := -1/240*x + 1/40;
#tau1 := -1/240*x^3 + 11/40 + 9/80*x;
#checkWeifeng(sigma1, [g1*g2*g3*g4], x, "PASSED");
#checkWeifeng(tau1, [g1*g2*g3*g4], x, "PASSED");
#sigma3 := 1/12*x^2-1/4*x+1/6+2113829/1318494300*(x+2)*(x+1)*(-x+3);
#tau3 := 51766519/659247150-2113829/1318494300*x;
#checkWeifeng(sigma3, [g1*g2*g4], x, "PASSED");
#checkWeifeng(tau3, [g1*g2*g4], x, "PASSED");

#checkWeifeng(-1, [-(x+1), x-1], x, "PASSED");
#checkWeifeng((x+1)*(x-1), [-1, x, x^2], x, "PASSED");

checkWeifeng(x+2, [(x+1), -(x-1)], x, "PASSED");
