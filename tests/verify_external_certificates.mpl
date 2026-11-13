# Run from the root of a source tree with the same flags as the benchmark.
# Recomputes spCertificates for each case of the external suites named in
# the BENCH_SUITES environment variable (comma separated, e.g. "left,arch")
# and checks the certificates exactly. VERIFY_TIMEOUT (default 300) limits
# the certificate computation in seconds. Prints one line per case:
#   VERIFY <suite> <case> <identity ok> <coefficients nonnegative> <degrees>
# or VERIFY <suite> <case> timeout|error|check_timeout.
interface(echo=0, warnlevel=0):
$include "src/StrictlyPositiveCert.mpl";
limit_text := cat("", getenv("VERIFY_TIMEOUT")):
compute_limit := `if`(limit_text = "", 300, parse(limit_text)):
for suite in StringTools:-Split(getenv("BENCH_SUITES"), ",") do
    read cat("tests/benchmarks/", suite, ".mpl");
    for k from 1 to nops(external_benchmarks) do
        b := external_benchmarks[k];
        try
            cert := timelimit(compute_limit, StrictlyPositiveCert:-spCertificates(b[1], b[2], x));
        catch "time expired":
            printf("VERIFY %s %d timeout\n", suite, k);
            next;
        catch:
            printf("VERIFY %s %d error\n", suite, k);
            next;
        end try;
        try
            printf("VERIFY %s %d %a %a %a\n", suite, k,
                timelimit(600, evalb(expand(b[1] - StrictlyPositiveCert:-dot_product([1, op(b[2])], cert)) = 0)),
                timelimit(600, andmap(h -> SolveTools:-SemiAlgebraic([h < 0], [x]) = [], cert)),
                map(h -> degree(expand(h), x), cert));
        catch:
            printf("VERIFY %s %d check_timeout %a\n", suite, k, map(h -> degree(expand(h), x), cert));
        end try;
    end do;
end do:
quit:
