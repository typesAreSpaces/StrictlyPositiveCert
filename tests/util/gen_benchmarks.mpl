_pwd := currentdir();
currentdir(FileTools[ParentDirectory](__FILE__));

read "nat_gens.mpl";

#read "gen_benchmark.mpl";
read "gen_benchmark_2.mpl";
read "random_poly_gen.mpl";

currentdir(_pwd);

leftStrictPositiveBenchmark := proc(num_iterations, basisGeneratorFunc)
local i;
local intervals1, basis1, isolated_points1, archimedean_poly1;
local f1;
local benchmarks := [];

    for i from 1 to num_iterations do
# -----------------------------------------------------------------------------------------------------
        intervals1 := intervalsGenerator(x, 3, 1, 10, 100);
        basis1 := basisGeneratorFunc(x, intervals1);
        isolated_points1 := getIsolatedPoints(intervals1);
        archimedean_poly1 := getArchimedeanPolynomial(x, isolated_points1, 2);

        f1 := strictlyLeftPolynomial(x, isolated_points1, 10);
        benchmarks := [op(benchmarks), [f1, basis1, archimedean_poly1]];
# -----------------------------------------------------------------------------------------------------
    end do;
    return benchmarks;
end proc;

rightStrictPositiveBenchmark := proc(num_iterations, basisGeneratorFunc)
local i;
local test_name1;
local intervals1, basis1, isolated_points1, archimedean_poly1;
local f1;
local benchmarks := [];

    for i from 1 to num_iterations do
# -----------------------------------------------------------------------------------------------------
        intervals1 := intervalsGenerator(x, 3, 1, 10, 100);
        basis1 := basisGeneratorFunc(x, intervals1);
        isolated_points1 := getIsolatedPoints(intervals1);
        archimedean_poly1 := getArchimedeanPolynomial(x, isolated_points1, 2);

        f1 := strictlyRightPolynomial(x, isolated_points1, 10);
        benchmarks := [op(benchmarks), [f1, basis1, archimedean_poly1]];
# -----------------------------------------------------------------------------------------------------
    end do;
    return benchmarks;
end proc;

productGensBenchmark := proc(num_iterations, basisGeneratorFunc)
local i;
local test_name1;
local intervals1, basis1, isolated_points1, archimedean_poly1;
local f1;
local random_point;
local benchmarks := [];

    for i from 1 to num_iterations do
# -----------------------------------------------------------------------------------------------------
        #intervals1 := intervalsGenerator(x, 1, 1, 10, 100);
        random_point := genRandomPoint(1, 10, 100);
        #basis1 := basisGeneratorFunc(x, intervals1);
        basis1 := [x-random_point, -(x-random_point)];
        #isolated_points1 := getIsolatedPoints(intervals1);
        isolated_points1 := [random_point];
        archimedean_poly1 := getArchimedeanPolynomial(x, isolated_points1, 2);

        f1 := basis1[1]*basis1[2] + 35;
        benchmarks := [op(benchmarks), [f1, basis1, archimedean_poly1]];
# -----------------------------------------------------------------------------------------------------
    end do;
    return benchmarks;
end proc;

archimedeanBenchmark:= proc(num_iterations, basisGeneratorFunc)
local i;
local test_name1;
local intervals1, basis1, isolated_points1, archimedean_poly1;
local f1;
local benchmarks := [];

    for i from 1 to num_iterations do
# -----------------------------------------------------------------------------------------------------
        intervals1 := intervalsGenerator(x, 3, 1, 10, 100);
        #intervals1 := intervalsGenerator(x, 3, 1, 10, 50);
        basis1 := basisGeneratorFunc(x, intervals1);
        isolated_points1 := getIsolatedPoints(intervals1);
        archimedean_poly1 := getArchimedeanPolynomial(x, isolated_points1, 2);

        f1 := archimedean_poly1 + 10;
        benchmarks := [op(benchmarks), [f1, basis1, archimedean_poly1]];
# -----------------------------------------------------------------------------------------------------
    end do;
    return benchmarks;
end proc;
