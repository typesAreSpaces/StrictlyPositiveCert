$define AVERKOV_SOS(g, N)         ((g - gamma)/(gamma + eps))^(2*N)
$define AVERKOV_1_EXPR(g, N)      AVERKOV_SOS(g, N)*g
$define AVERKOV_EXPR(N)           add(term, term in map(g_i -> AVERKOV_1_EXPR(g_i, N), gens))

$ifdef ENABLE_DEBUGGING
$define DEBUG(F, L, X) lprint(">> Debugging file ", F, " at line ", L); X
$else
$define DEBUG(F, L, X)
$endif

# Timed routines pair INIT_START_LOG_TIME with END_LOG_TIME in a conditional
# finally block, covering returns and catchable errors. An outer timelimit or
# user interrupt can bypass finally blocks inside the interrupted computation.
$define START_LOG_TIME(X, S) stack_level:=stack_level+1;fd := FileTools:-Text:-Open("log_time.txt", append);local _log_time_S := time();FileTools:-Text:-WriteString(fd, cat("Start: ", X, " ", convert(stack_level, string), "\n"));FileTools:-Text:-Close(fd);
$define END_LOG_TIME(X, S) fd := FileTools:-Text:-Open("log_time.txt", append);FileTools:-Text:-WriteString(fd, cat("End: ", X, " ", convert(stack_level, string), "\nTime: ", convert(time() - _log_time_S, string), "\n"));FileTools:-Text:-Close(fd);stack_level:=stack_level-1;
$define INIT_START_LOG_TIME(X, S) local fd;START_LOG_TIME(X, S)

with(SolveTools, SemiAlgebraic);
with(RootFinding, Isolate);

StrictlyPositiveCert := module() option package;

export isSOS;
export quickCertificates;
export dot_product;
export bound_info;
export spCertificates;
export monogenicCertificates;
export SplitNonnegativeFactor;

$ifdef LOG_TIME
local stack_level := -1;
$endif

$include "src/utilities.mpl";
$include "src/univariate_signs.mpl";
$include "src/signature_constructions.mpl";
$include "src/bound_poly.mpl";
$include "src/eps_computation.mpl";
$include "src/averkov_constructions.mpl";
$include "src/monogenic_certificates.mpl";
$include "src/spCertificates.mpl";

end module;
