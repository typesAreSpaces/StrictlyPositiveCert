# Run from the repository root with the current Makefile's flags and values.
# Load source directly so the installed .mla does not affect the comparison.
interface(echo=0):
$include "src/StrictlyPositiveCert.mpl";
currentdir("tests"):
try
    read "benchmarking.mpl";
    printf("FACTOR_SOS_BENCHMARKING_COMPLETE\n");
catch:
    lprint("FACTOR_SOS_BENCHMARKING_ERROR", lastexception);
end try:
quit:
