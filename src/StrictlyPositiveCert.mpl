$define ENABLE_DEBUGGING          true
$define ENABLE_VERIFICATION       false
$define ENABLE_BINARY_SEARCH_AVKL true
$define ENABLE_BINARY_SEARCH_LS   false
$define ENABLE_N_HEURISTIC        false
$define ENABLE_AVERKOV_CHECK      false
#$define EPS_FACTOR               17/10
#$define EPS_FACTOR               1/10
$define EPS_FACTOR                1/100
$define ENABLE_POST_OPT           false
#$define WEIFENG_OPTIMIZATION
$define N_GUESS_AVKL              70
$define N_GUESS_LS                600
#$define LOG_TIME
$define SOS_DELAY_SEARCH          5

$define DEBUG_EXIT lprint(">> Debugging, getting out"); return 0
$define DEBUG(F, L, y, x) if (y) then lprint(">> Debugging file ", F, " at line ", L); x; end if

$define START_LOG_TIME(X, S) stack_level:=stack_level+1;fd := FileTools:-Text:-Open("log_time.txt", append);local _log_time_S := time();FileTools:-Text:-WriteString(fd, cat("Start: ", X, " ", convert(stack_level, string), "\n"));FileTools:-Text:-Close(fd);
$define END_LOG_TIME(X, S) fd := FileTools:-Text:-Open("log_time.txt", append);FileTools:-Text:-WriteString(fd, cat("End: ", X, " ", convert(stack_level, string), "\nTime: ", convert(time() - _log_time_S, string), "\n"));FileTools:-Text:-Close(fd);stack_level:=stack_level-1;
$define INIT_START_LOG_TIME(X, S) local fd;START_LOG_TIME(X, S)

with(SolveTools, SemiAlgebraic);
with(RootFinding, Isolate);
with(Optimization, Minimize);
with(RegularChains, SemiAlgebraicSetTools, PolynomialRing);

StrictlyPositiveCert := module() option package;

export dot_product;
export bound_info;
export spCertificates;

$ifdef LOG_TIME
local stack_level := -1;
$endif

$include "src/utilities.mpl";
$include "src/signature_constructions.mpl";
$include "src/bound_poly.mpl";
$include "src/lower_bound_poly.mpl";
$include "src/eps_computation.mpl";
$include "src/averkov_constructions.mpl";
$include "src/spCertificates.mpl";

end module;
