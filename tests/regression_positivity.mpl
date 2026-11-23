# Read current source rather than a compiled archive.
read "src/StrictlyPositiveCert.mpl":
kernelopts(opaquemodules=false):

check_positivity_traversal := proc()
local x, item, answer, cells, interval, bounds;
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
        # Both endpoint conventions, repeated roots, and an interior sign dip.
        [[[0<=x,x<=1]], x, false],
        [[[0<=x,x<=1]], 1-x, false],
        [[[1<=x,x<=2]], (x^2-2)^2, false],
        [[[1<=x,x<=2]], (x^2-2)^2-1/10^20, false],
        [[[1<=x,x<=2]], (x^2-2)^2+1/10^20, true],
        # Algebraic endpoints and isolated points use exact signs.
        [[[sqrt(2)<=x,x<=2]], x^2-2, false],
        [[[sqrt(2)<=x,x<=2]], x^2-1, true],
        [[[x=sqrt(2)]], x^2-1, true],
        [[[x=sqrt(2)]], x^2-2, false],
        [[[x=0]], 0, false],
        [[[x=0]], -1, false],
        [[[x=0]], 1, true],
        # Rays and the whole line must not substitute infinity into p.
        [[[x<=-1]], -x, true],
        [[[x<=1]], x^2-2, false],
        [[[1<=x]], x, true],
        [[[x=x]], x^2+1, true],
        [[[x=x]], (x-1)^2, false],
        [[], -1, true]
    ] do
        # Convert legacy interval fixtures here; the production pipeline uses cells.
        cells := [];
        for interval in item[1] do
            bounds := StrictlyPositiveCert:-bound_info(x,interval,0);
            cells := [op(cells),[op(bounds),evalb(bounds[1]<>-infinity),
                                          evalb(bounds[2]<>infinity)]];
        end do;
        answer := StrictlyPositiveCert:-univariate_positive_on_cells(item[2],cells,x);
        if answer <> item[3] then
            error "incorrect interval positivity result for %1: %2", item, answer;
        end if;
    end do;
    printf("PASS: positivity interval and critical-point traversal\n");
end proc:
check_positivity_traversal():
