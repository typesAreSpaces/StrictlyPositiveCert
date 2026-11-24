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
$define INIT_START_LOG_TIME(X) local fd;stack_level:=stack_level+1;fd := FileTools:-Text:-Open("log_time.txt", append);local _log_time_S := time();log_entries := [op(log_entries), [X, _log_time_S]];FileTools:-Text:-WriteString(fd, cat("Start: ", X, " ", convert(stack_level, string), "\n"));FileTools:-Text:-Close(fd);
$define END_LOG_TIME(X) fd := FileTools:-Text:-Open("log_time.txt", append);FileTools:-Text:-WriteString(fd, cat("End: ", X, " ", convert(stack_level, string), "\nTime: ", convert(time() - _log_time_S, string), "\n"));FileTools:-Text:-Close(fd);log_entries := log_entries[1..-2];stack_level:=stack_level-1;

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
# Open timing entries [routine, start time], innermost last.
local log_entries := [];

# Close entries a timelimit abandoned: Maple skips finally blocks below the
# timelimit boundary, so a caller above it writes their End records.
local log_unwind := proc(level)
local fd, entry;
    while stack_level > level and log_entries <> [] do
        entry := log_entries[-1];
        fd := FileTools:-Text:-Open("log_time.txt", append);
        FileTools:-Text:-WriteString(fd, cat("End: ", entry[1], " ", convert(stack_level, string),
                                             "\nTime: ", convert(time() - entry[2], string), "\n"));
        FileTools:-Text:-Close(fd);
        log_entries := log_entries[1..-2];
        stack_level := stack_level-1;
    end do;
end proc;
$endif

$include "src/utilities.mpl";
$include "src/univariate_signs.mpl";
$include "src/signature_constructions.mpl";
$include "src/bound_poly.mpl";
$include "src/eps_computation.mpl";
$include "src/univariate_heuristics.mpl";
$include "src/averkov_constructions.mpl";
$include "src/monogenic_certificates.mpl";
$include "src/spCertificates.mpl";

end module;
