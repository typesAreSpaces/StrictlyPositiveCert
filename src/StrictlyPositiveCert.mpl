# Multiplier of averkov_lemma_7 for a generator g at exponent N. The default
# is T_N((g-gamma)/gamma)^2/T_N(1+eps/gamma)^2 (Chebyshev); -DAVERKOV_POWER restores
# Averkov's ((g-gamma)/(gamma+eps))^(2N). gamma, eps, lambda are locals there.
$ifdef AVERKOV_POWER
$define AVERKOV_SOS(g, N)         ((g - gamma)/(gamma + eps))^(2*N)
$else
$define AVERKOV_SOS(g, N)         (lambda(N)*averkov_chebyshev(N, (g - gamma)/gamma)^2)
$endif
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
$define START_LOG_TIME(X) stack_level:=stack_level+1;fd := FileTools:-Text:-Open("log_time.txt", append);local _log_time_start := time();FileTools:-Text:-WriteString(fd, cat("Start: ", X, " ", convert(stack_level, string), "\n"));FileTools:-Text:-Close(fd);
$define END_LOG_TIME(X) fd := FileTools:-Text:-Open("log_time.txt", append);FileTools:-Text:-WriteString(fd, cat("End: ", X, " ", convert(stack_level, string), "\nTime: ", convert(time() - _log_time_start, string), "\n"));FileTools:-Text:-Close(fd);stack_level:=stack_level-1;
$define INIT_START_LOG_TIME(X) local fd;START_LOG_TIME(X)

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

# Saved in debug builds so the message identifies the library's source branch.
$ifdef ENABLE_DEBUGGING
local ModuleLoad := proc()
    printf("StrictlyPositiveCert loaded (branch: fixing_lower_bound_poly_no_ext_lemma_eps_no_avkl)\n");
end proc;
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
