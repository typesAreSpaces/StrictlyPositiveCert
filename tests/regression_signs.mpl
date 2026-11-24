# Compare the specialized exact routines with the general solver as an oracle.
read "src/StrictlyPositiveCert.mpl":
kernelopts(opaquemodules=false):

check_sign_routines := proc()
local x, p, q, polynomials, expected, actual, S, cells, cell, strict, condition;
local chart, relation, oracle, result, i, item, h, g, big, small, status, bounds, pieces;
local m, M, queries := 0;
    polynomials := [0,1,-1,x,(x-1)^2,-(x-1)^2,x^2-2,
                    (x^2-2)^2-1/100,x^5-x-1,x*(x^2-2),
                    (100000*x-1)*(100000*x-2)];
    for p in polynomials do
        actual := StrictlyPositiveCert:-univariate_nonnegative_on_R(p,x);
        expected := evalb(SolveTools:-SemiAlgebraic([p<0],[x])=[]);
        if actual <> expected then error "global nonnegativity mismatch: %1",p; end if;
        for q in [1,-1,1-x^2,-x^2,x^2-2] do
            S := StrictlyPositiveCert:-univariate_sign_set([p,q],
                                                           ["nonnegative","nonnegative"],x);
            oracle := SolveTools:-SemiAlgebraic([p>=0,q>=0],[x]);
            actual := map(piece -> StrictlyPositiveCert:-bound_info(x,piece,0),S);
            expected := map(piece -> StrictlyPositiveCert:-bound_info(x,piece,0),oracle);
            expected := sort(expected,(a,b)->StrictlyPositiveCert:-univariate_compare(a[1],b[1])<0);
            if nops(actual)<>nops(expected) then
                error "component count mismatch: %1",[p,q,S,oracle];
            end if;
            for i from 1 to nops(actual) do
                if StrictlyPositiveCert:-univariate_compare(actual[i][1],expected[i][1])<>0 or
                StrictlyPositiveCert:-univariate_compare(actual[i][2],expected[i][2])<>0 then
                    error "algebraic boundary mismatch: %1",[p,q,S,oracle];
                end if;
            end do;
            queries := queries+1;
        end do;
        # Include open/closed endpoints, singleton algebraic points, and rays.
        for cells in [[],[[-1,1,true,true]],[[-1,1,false,false]],
                      [[1,2,true,false]],[[sqrt(2),sqrt(2),true,true]],
                      [[-infinity,0,false,false],[1,infinity,true,false]]] do
            for strict in [true,false] do
                actual := StrictlyPositiveCert:-univariate_positive_on_cells(p,cells,x,strict);
                expected := true;
                condition := `if`(strict,p<=0,p<0);
                if cells=[[sqrt(2),sqrt(2),true,true]] then
                    # SemiAlgebraic accepts rational input coefficients only.
                    pieces := [[x^2=2,x>0]];
                else
                    pieces := StrictlyPositiveCert:-univariate_cells_to_sas(cells,x);
                end if;
                for cell in pieces do
                    if SolveTools:-SemiAlgebraic([op(cell),condition],[x])<>[] then expected:=false; break; end if;
                end do;
                if actual<>expected then
                    error "cell positivity mismatch: %1",[p,cells,strict,actual,expected];
                end if;
                queries := queries+1;
            end do;
        end do;
    end do;
    for g in [1-x^2,-x^2,-1,-(x^2-1)*(x^2-4)] do
        for h in [0,1,-1,x,x^2,(x-1)^2,x^2-2,(x^2-2)^2] do
            if SolveTools:-SemiAlgebraic([h<=0],[x])=[] then status:="ok";
            else
                big:=evalb(SolveTools:-SemiAlgebraic([g>0,h<=0],[x])<>[]);
                small:=evalb(SolveTools:-SemiAlgebraic([g<=0,h<=0],[x])<>[]);
                if big and small then status:="infeasible";
                elif big then status:="big"; else status:="small"; end if;
            end if;
            actual := StrictlyPositiveCert:-monogenic_status_lambda(h,g,0,1,x);
            if actual<>status then error "lambda classification mismatch: %1",[h,g,actual,status]; end if;
            queries:=queries+1;
        end do;
    end do;
    # Preserve the original dyadic m/M values, including equality and M's cap.
    for item in [[x+3/2,1-x^2],[x+1,-x^2],[0,-1],
                 [1/8,1/1000-x^2],[32,8-x^2],[x^2-1/2,-(x^2-1)*(x^2-4)]] do
        p,g := op(item); m:=1; M:=1;
        if SolveTools:-SemiAlgebraic([g>=0],[x])<>[] then
            while SolveTools:-SemiAlgebraic([g>=0,p<2*m],[x])=[] do m:=2*m; end do;
            while SolveTools:-SemiAlgebraic([g>=0,p<m],[x])<>[] do m:=m/2; end do;
            while SolveTools:-SemiAlgebraic([g>=0,g>M],[x])<>[] do M:=2*M; end do;
            while M>1/2^20 and SolveTools:-SemiAlgebraic([g>=0,g>M/2],[x])=[] do M:=M/2; end do;
        end if;
        bounds:=StrictlyPositiveCert:-monogenic_search_bounds(p,g,x);
        if bounds[1..2]<>[m,M] then error "normalization changed: %1",[item,bounds,[m,M]]; end if;
    end do;
    # Root isolation against algebraic cell endpoints: roots at, near, and
    # beyond sqrt(2) and sqrt(3), including a root 10^-17 from sqrt(2) that
    # no bounded refinement separates and must reach the Sturm fallback.
    for g in [2-x^2,-(x^2-2)*(x^2-3)] do
        cells := StrictlyPositiveCert:-univariate_select_cells(
                     StrictlyPositiveCert:-univariate_sign_chart(g,x),"nonnegative");
        for p in [x^2-2,(x^2-2)^2,2-x^2,(x^2-2)*(x-1),(x^2-2)^2*(x-3/2),
                  x-14142/10000,x-14143/10000,(x-14142/10000)*(x-14143/10000),
                  x-14142135623730951/10^16,x-14142135623730950/10^16,
                  (x^2-3)*(x+2),(x^2-2)*(x^2-3)+1/10^12,x^2+1] do
            for strict in [true,false] do
                actual := StrictlyPositiveCert:-univariate_positive_on_cells(p,cells,x,strict);
                expected := evalb(SolveTools:-SemiAlgebraic([g>=0,`if`(strict,p<=0,p<0)],[x])=[]);
                if actual<>expected then
                    error "algebraic endpoint positivity mismatch: %1",[p,g,strict,actual,expected];
                end if;
                queries := queries+1;
            end do;
        end do;
    end do;
    cells := StrictlyPositiveCert:-univariate_select_cells(
                 StrictlyPositiveCert:-univariate_sign_chart(2-x^2,x),"nonnegative");
    if StrictlyPositiveCert:-univariate_roots_in_interiors(
           10^16*x-14142135623730951,cells,x,4) <> FAIL then
        error "near-endpoint root was decided without enough refinement";
    end if;
    printf("PASS: %d sign/set/classification oracle comparisons and exact dyadic bounds\n",queries);
end proc:
check_sign_routines():
