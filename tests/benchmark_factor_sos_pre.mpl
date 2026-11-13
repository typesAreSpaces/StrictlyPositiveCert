# Run from the repository root with the same flag combinations as the tests.
# Reports warm, per-call wall times; compilation and validation are excluded.
$define FACTOR_SOS_BENCHMARK
$include "tests/factor_sos_pre.mpl";

benchmark := proc(label, operation, certificate := false)
local result, reps := 1, samples := [], elapsed, started, batch, i, degrees := "-";
local factor_flag := 0, weifeng_flag := 0;
$ifdef FACTOR_SOS_PRE
    factor_flag := 1;
$endif
$ifdef WEIFENG_OPTIMIZATION
    weifeng_flag := 1;
$endif
    result := operation();
    if certificate then
        degrees := StringTools:-Join(map(d -> convert(d, string),
            map(h -> degree(expand(h), x), result)), ";");
    end if;
    # Calibrate batches to at least 50 ms to avoid timer resolution artifacts.
    do
        started := time[real]();
        for i to reps do result := operation(); end do;
        elapsed := time[real]()-started;
        if elapsed >= 0.05 or reps >= 4096 then break; end if;
        reps := 2*reps;
    end do;
    for batch to 7 do
        started := time[real]();
        for i to reps do result := operation(); end do;
        samples := [op(samples), 1000*(time[real]()-started)/reps];
    end do;
    samples := sort(samples);
    printf("BENCH,%s,%d,%d,%d,%.6f,%.6f,%.6f,%s\n",
        label, factor_flag, weifeng_flag, reps,
        samples[4], samples[1], samples[7], degrees);
end proc:

runBenchmarks := proc()
local S, basis, f, test_case;
    runTests();
    benchmark("regression_suite", () -> runTests());
    benchmark("isSOS_positive", () -> StrictlyPositiveCert:-isSOS((x^2+1)^3, x));
    benchmark("isSOS_sign_change", () -> StrictlyPositiveCert:-isSOS((x^2+1)*(x^2-1), x));
    S := SemiAlgebraic([x >= -2, x <= 2], [x]);
    benchmark("interval_positive", () -> StrictlyPositiveCert:-checkPositivityOverSAS(S, (x^2+1)^3, x));
    benchmark("interval_sign_change", () -> StrictlyPositiveCert:-checkPositivityOverSAS(S, (x^2+1)*(x^2-1), x));
    benchmark("minimum", () -> StrictlyPositiveCert:-computeMin(S, (x^2+1)^3+7, x));
    benchmark("epsilon", () -> StrictlyPositiveCert:-findEps(x,
        [1-x^2, (x^2+2)*(1-x^2)], (x^2+1)*(2-x^2)));
    benchmark("lower_bound_constant", () -> StrictlyPositiveCert:-find_constant_lower_bound_poly(
        (x^2+1)^2, 1, 1-x^2, x));
    benchmark("lower_bound", () -> [StrictlyPositiveCert:-lower_bound_poly(x,
        (x^2+1)^2*(x+2), 1-x^2)]);

    basis := [-x^2/2, -x^2/2];
$ifdef WEIFENG_OPTIMIZATION
    f := (x^2+1)*(-x^2/2+(x^2+2)*(x^2+3)*(1-x^2));
$else
    f := (x^2+1)*(-x^2/2+(x^2+2)*(-x^4+(x^2+3)*(1-x^2)));
$endif
    benchmark("controlled_stages", () -> Stages:-spCertificates(f, basis, x), true);
    for test_case in [
        ["certificate_factored_input", (x^2+1)*(x+2), [1-x^2]],
        ["certificate_linear", x+1, [-x^2]],
        ["certificate_global_SOS", (x^2+1)*(x^2+2), [-x^2]],
        ["certificate_minus_one", (x^2+1)*(x+1), [-1]]
    ] do
        f := test_case[2]; basis := test_case[3];
        benchmark(test_case[1], () -> StrictlyPositiveCert:-spCertificates(f, basis, x), true);
    end do;
    printf("BENCHMARK_COMPLETE\n");
end proc:
runBenchmarks():
quit:
