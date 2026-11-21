# Run from the repository root with the same -D values as the Makefile.
# Read source directly so this does not depend on a rebuilt .mla library.
read "src/StrictlyPositiveCert.mpl":
kernelopts(opaquemodules=false):

check_epsilon := proc()
local x;
    # Shared negative intersection roots must be counted after deduplication.
    # On {f <= 0} = (-infinity, 3/4], max(min(basis)) = -1/4.
    if StrictlyPositiveCert:-findEps(x,[x-1,2*x-1,3*x-1,4-x^2],x-3/4) <> 1/4 then
        error "epsilon computation failed for coincident generator intersections";
    end if;
    printf("PASS: duplicate-intersection epsilon\n");
end proc:
check_epsilon():
