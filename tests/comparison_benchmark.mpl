with(StrictlyPositiveCert, spCertificates);

interface(echo=0);
getenv(OMP_NUM_THREADS);
kernelopts(numcpus);

printlevel := -1;
#printlevel := 8;

read "util/check_weifeng.mpl";
read "util/check_realcertify.mpl";

# Load RealCertify
currentdir(homedir);
currentdir("Documents/GithubProjects/RealCertify");
read "multivsos/multivsos.mm";

left_comparison_benchmark := proc(n)
local i, x, f, g;

    f := x + n;
    g := simplify(-product((x^2 - i^2)^2, i = 1 .. n)/(x + n)^2);

    #try
        #checkRealCertify((x+n)^2*f, [(x+n)^2*g], "Comparison benchmark using RealCertify");
    #catch:
        #printf(">> RealCertify fails\n");
    #end try;
    #try
        #checkRealCertify(f, [g], "Comparison benchmark using RealCertify");
    #catch:
        #printf(">> RealCertify fails\n");
    #end try;
    try
        checkWeifeng(f, [g], x, "Comparison benchmark using approach");
    catch:
        printf(">> Our approach fails\n");
    end try;
end proc;

in_between_comparison_benchmark := proc(n)
local i, x, f, g;

    f := (x+1)*(x-1);
    g := simplify(-product((x^2 - i^2)^2, i = 1 .. n)/(x + 1)^2/(x - 1)^2);

    #try
        #checkRealCertify((x+n)^2*f, [(x+n)^2*g], "Comparison benchmark using RealCertify");
    #catch:
        #printf(">> RealCertify fails\n");
    #end try;
    #try
        #checkRealCertify(f, [g], "Comparison benchmark using RealCertify");
    #catch:
        #printf(">> RealCertify fails\n");
    #end try;
    try
        checkWeifeng(f, [g], x, "Comparison benchmark using approach");
    catch:
        printf(">> Our approach fails\n");
    end try;
end proc;

# --------------------
# Comparison Benchmark
# --------------------

left_comparison_benchmark(1);

printf("\n>> Start benchmark\n");

#checkWeifeng(-(x-1), [-(x+2)*x], x, "PASSED");

local n_max := 9;
for i from 1 to n_max do
  in_between_comparison_benchmark(i);
end do;
