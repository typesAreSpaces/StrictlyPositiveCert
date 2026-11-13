# Run from the root of a source tree. Times spCertificates (via checkWeifeng)
# on every case of tests/benchmarks/lbp_stress.mpl, with a 60 s time limit per
# case. Cases listed in the SKIP_CASES environment variable (comma separated)
# are reported as skipped without running them.
interface(echo=0, warnlevel=0):
$include "src/StrictlyPositiveCert.mpl";
currentdir("tests"):
read "util/check_weifeng.mpl":
read "benchmarks/lbp_stress.mpl":
skip_text := cat("", getenv("SKIP_CASES")):
skip := `if`(skip_text = "", {}, {op(map(parse, StringTools:-Split(skip_text, ",")))}):
for k from 1 to nops(external_benchmarks) do
    b := external_benchmarks[k];
    if member(k, skip) then
        printf("\n>> Test\nBatch: %d Test - Weifeng\n>> Skipped after timeout\n", k);
        next;
    end if;
    try
        timelimit(60, checkWeifeng(b[1], b[2], x, cat("Batch: ", k, " Test - Weifeng")));
    catch "time expired":
        printf(">> Timeout Weifeng\n");
    catch:
        printf(">> Error Weifeng %s\n", StringTools:-FormatMessage(lastexception[2..-1]));
    end try;
end do:
printf("\nLBP_STRESS_COMPLETE\n"):
quit:
