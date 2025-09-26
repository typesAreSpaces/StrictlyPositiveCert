with(StrictlyPositiveCert, spCertificates);

benchmarking_pwd := currentdir();

interface(echo=0);
getenv(OMP_NUM_THREADS);
kernelopts(numcpus);

printlevel := -1;
#printlevel := 8;

read "util/gen_benchmarks.mpl";
read "util/check_weifeng.mpl";
read "util/check_realcertify.mpl";

# Load RealCertify
currentdir(homedir);
currentdir("Documents/GithubProjects/RealCertify");
read "multivsos/multivsos.mm";interface(echo=0);

natGensBenchmark := proc(num_iterations, basisGeneratorFunc, max_time)
local i;
local test_name1, test_name2, test_name3, test_name4;
local intervals1, basis1, isolated_points1;
local nat_gens;

    for i from 1 to num_iterations do
# -----------------------------------------------------------------------------------------------------
        intervals1 := intervalsGenerator(x, 3, 1, 10, 100);
        basis1 := basisGeneratorFunc(x, intervals1);
        print("Basis", basis1);
        isolated_points1 := getIsolatedPoints(intervals1);
        nat_gens := gen_nat_gens(basis1, x);
# -----------------------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------------------
# Testing Weifeng certificates
        test_name1:=cat("Batch: ", convert(i, string),  " Test Weifeng");
        try
            map(f1 -> timelimit(max_time, checkWeifeng(f1, basis1, x, test_name1)), nat_gens);
            printf(">> Succeeds Weifeng %s\n", test_name1);
        catch:
            printf(">> Timeout Weifeng %s\n", test_name1);
        end try;
# -----------------------------------------------------------------------------------------------------
    end do;
end proc;

natGensBenchmark2 := proc(x, num_points)
local i, _p;

local pos_isolated_points := [seq(i, i=1 .. num_points)];
local isolated_points := sort([op(map(p -> -p, pos_isolated_points)), op(pos_isolated_points)]);
local basis_poly := mul(_p, _p in map(p -> x - p, isolated_points));
local basis := [-basis_poly, basis_poly];
local tempPoly;
    printf("\n>> Start benchmark\n");

# Left natural generator
    try
        checkWeifeng(x + num_points, basis, x, "Left Natural Generator Test");
    catch:
        printf(">> Timeout Weifeng %s\n", "Left Natural Generator Test");
    end try;
# In between natural generators
    for i from 1 to 2*num_points - 1 do
        tempPoly := (x - isolated_points[i])*(x - isolated_points[i+1]);
        try
            checkWeifeng(tempPoly, basis, x, cat("In between Natural Generator Test ", convert(i, string)));
        catch:
            printf(">> Timeout Weifeng %s\n", cat("In between Natural Generator Test ", convert(i, string)));
        end try;
    end do;
# Right natural generator
    try
        checkWeifeng(-(x - num_points), basis, x, "Right Natural Generator Test");
    catch:
        printf(">> Timeout Weifeng %s\n", "Right Natural Generator Test");
    end try;

# -----------------------------------------------------------------------------------------------------
end proc;

runBenchmark := proc(benchmarks, max_time)
    printf("\n>> Start benchmark\n");
local i := 1;
local test_name;
    for benchmark in benchmarks do
# -----------------------------------------------------------------------------------------------------
# Testing RealCertify
        test_name:=cat("Batch: ", convert(i, string),  " Test - RealCertify without Archimedean polynomial");

        try
            timelimit(max_time, checkRealCertify(benchmark[1], benchmark[2], test_name));
            printf(">> Succeeds RealCertify %s\n", test_name);
        catch:
            printf(">> Timeout RealCertify %s\n", test_name);
        end try;

        test_name:=cat("Batch: ", convert(i, string),  " Test - RealCertify with Archimedean polynomial");

        try
            timelimit(max_time, checkRealCertify(benchmark[1], [op(benchmark[2]), benchmark[3]], test_name));
            printf(">> Succeeds RealCertify %s\n", test_name);
        catch:
            printf(">> Timeout RealCertify %s\n", test_name);
        end try;
# -----------------------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------------------
# Testing Weifeng certificates
        test_name:=cat("Batch: ", convert(i, string),  " Test - Weifeng");

        try
            timelimit(max_time, checkWeifeng(benchmark[1], benchmark[2], x, test_name));
            printf(">> Succeeds Weifeng %s\n", test_name);
        catch:
            printf(">> Timeout Weifeng %s\n", test_name);
        end try;
# -----------------------------------------------------------------------------------------------------
        i := i + 1;
    end do;
end proc;

# This assumes isolated points case
runBenchmark1 := proc()
#Status PASSED
#Heuristics: -D USE_BOUNDS
    for d from 1 to 7 do:
      natGensBenchmark2(x, d);
    end do; 
end proc;

certificates_of_strictly_positive_polynomials_benchmark := proc()
local num_iterations := 100;
local max_time := 300;
#Status PASSED
#Heuristics: No heuristics, -D USE_BOUNDS
    runBenchmark(leftStrictPositiveBenchmark(num_iterations, linearBasisNoBoundedPoly), max_time);

#Status PASSED
#Heuristics: No heuristics, -D USE_BOUNDS
    runBenchmark(rightStrictPositiveBenchmark(num_iterations, linearBasisNoBoundedPoly), max_time);

#Status PASSED
#Heuristics: No heuristics, -D USE_BOUNDS
    runBenchmark(productGensBenchmark(num_iterations, linearBasisNoBoundedPoly), max_time);

#Status PASSED
#Heuristics: -D USE_BOUNDS
    runBenchmark(archimedeanBenchmark(num_iterations, linearBasisNoBoundedPoly), max_time);
end proc;

strictly_positive_polynomials_close_to_semialgebraic_sets_benchmark := proc()
  printf("\n>> Start benchmark\n");
local test_name := "Strictly positive polynomials close to semialgebraic sets";
local max_time := 300;

  for eps in [1/2, 1/3] do
    for k from 13 to 21 by 2 do
      local f := (1 + eps) + x;
      local g := [(1 - x^2)^k];

      try
        timelimit(max_time, checkRealCertify(f, g, test_name));
        printf(">> Succeeds RealCertify %s\n", test_name);
      catch:
        printf(">> Timeout RealCertify %s\n", test_name);
      end try;

      try
        timelimit(max_time, checkWeifeng(f, g, x, test_name));
        printf(">> Succeeds Weifeng %s\n", test_name);
      catch:
        printf(">> Timeout Weifeng %s\n", test_name);
      end try;
    end do;
  end do;
  return 0;
end proc;

strictly_positive_polynomials_stengle_benchmark := proc()
  printf("\n>> Start benchmark\n");
local test_name := "Strictly positive polynomials Stengles degree estimation";
local max_time := 300;

  for d in [10, 5, 1, 1/2, 1/3] do
    for k from 1 to 3 by 2 do
      local f := (1 - x^2) + d;
      local g := [(1 - x^2)^k];

      try
        timelimit(max_time, checkRealCertify(f, g, test_name));
        printf(">> Succeeds RealCertify %s\n", test_name);
      catch:
        printf(">> Timeout RealCertify %s\n", test_name);
      end try;

      try
        timelimit(max_time, checkWeifeng(f, g, x, test_name));
        printf(">> Succeeds Weifeng %s\n", test_name);
      catch:
        printf(">> Timeout Weifeng %s\n", test_name);
      end try;
    end do;
  end do;
  return 0;
end proc;


# Reads a file 'path' with a benchmark.
# Assumes the benchmark is stored in a variable called 'external_benchmarks'
runExternalBenchmark := proc(path)
  read path;
local max_time := 300;
    runBenchmark(external_benchmarks, max_time); 
end proc;

# ----------
# Benchmarks
# ----------

printf("\n>> Start Benchmarks\n");

#runBenchmark1();
#certificates_of_strictly_positive_polynomials_benchmark();
strictly_positive_polynomials_close_to_semialgebraic_sets_benchmark();
#strictly_positive_polynomials_stengle_benchmark();

#runExternalBenchmark(cat(benchmarking_pwd, "/benchmarks/left.mpl"));
#runExternalBenchmark(cat(benchmarking_pwd, "/benchmarks/right.mpl"));
#runExternalBenchmark(cat(benchmarking_pwd, "/benchmarks/lifted.mpl"));
#runExternalBenchmark(cat(benchmarking_pwd, "/benchmarks/arch.mpl"));

