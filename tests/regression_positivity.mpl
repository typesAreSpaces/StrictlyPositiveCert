# Read current source rather than a compiled archive.
read "src/StrictlyPositiveCert.mpl":
kernelopts(opaquemodules=false):

check_positivity_traversal := proc()
local x, item, answer;
    for item in [
        # Earlier derivative roots must not skip the minimum at x=2.
        [[[1<=x,x<=3]], (x+1)^2*(x-2)^2-1/100, false],
        # No derivative roots: every interval still needs endpoint checks.
        [[[-2<=x,x<=-1],[1<=x,x<=2]], -x, false],
        # The same applies after consuming the final critical point.
        [[[-1<=x,x<=1],[3<=x,x<=4]], 4-x^2, false],
        # A critical point beyond the right endpoint must not be evaluated.
        [[[-3<=x,x<=-2]], x^2-1, true],
        # Preserve critical points for later intervals, even in reverse order.
        [[[5<=x,x<=6],[1<=x,x<=3]], (x+1)^2*(x-2)^2-1/100, false],
        # Strict positivity rejects zeros, including at an interior minimum.
        [[[-1<=x,x<=1]], x^2, false],
        [[[-2<=x,x<=-1],[1<=x,x<=2]], x^2+1, true],
        [[], -1, true]
    ] do
        answer := StrictlyPositiveCert:-checkPositivityOverSAS(item[1],item[2],x);
        if answer <> item[3] then
            error "incorrect interval positivity result for %1: %2", item, answer;
        end if;
    end do;
    printf("PASS: positivity interval and critical-point traversal\n");
end proc:
check_positivity_traversal():
