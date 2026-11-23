read "src/StrictlyPositiveCert.mpl":
kernelopts(opaquemodules=false):
$ifdef LOG_TIME
if getenv("SPCERT_LOG_DIR")<>"" then currentdir(getenv("SPCERT_LOG_DIR")); end if:
$endif

check_univariate_heuristics := proc()
local x, B:=1-x^2, cells:=[[-1,1,true,true]], H, f, gens, weights, residual;
local result, caught, item;
local valid := proc(f,gens,H,x)
local h;
    if nops(H)<>nops(gens)+1 or expand(f-StrictlyPositiveCert:-dot_product([1,op(gens)],H))<>0 then
        error "certificate reconstruction failed: %1",[f,gens,H];
    end if;
    for h in H do
        if not type(h,polynom(rational,x)) or SolveTools:-SemiAlgebraic([h<0],[x])<>[] then
            error "invalid SOS multiplier: %1",h;
        end if;
    end do;
end proc;

    # A quotient with zeros is SOS, despite not being strictly positive.
    H:=StrictlyPositiveCert:-spCertificates(x^3+1/10,[B,x],x,max_averkov_N=0);
    if H<>[1/10,0,x^2] then error "division did not accept a square quotient"; end if;
    valid(x^3+1/10,[B,x],H,x);

    # A remainder positive only on B must continue to the monogenic stage.
    f:=x^5+x+3/2; gens:=[B,x^3];
    weights:=StrictlyPositiveCert:-univariate_reduction(f,gens,x,B,max_trials=0);
    if weights<>[0,x^2] then error "division reduction was not found"; end if;
    residual:=expand(f-StrictlyPositiveCert:-dot_product(gens,weights));
    if residual<>x+3/2 or SolveTools:-SemiAlgebraic([residual<0],[x])=[] then
        error "test remainder should not be globally nonnegative";
    end if;
    H:=StrictlyPositiveCert:-spCertificates(f,gens,x,max_averkov_N=0);
    valid(f,gens,H,x);

    # Cancellation between generators gives a constant direction. Coefficient
    # scaling reaches the useful rational lambda immediately, even for tiny gaps.
    gens:=[-x^2,x^2-1/10^6];
    H:=StrictlyPositiveCert:-spCertificates(-1,gens,x,max_averkov_N=0,max_heuristic_trials=4);
    if H<>[1,2*10^6,2*10^6] then error "constant direction did not use rational scaling"; end if;
    valid(-1,gens,H,x);
    H:=StrictlyPositiveCert:-spCertificates(-1,[1,gens[1],0,gens[2]],x,
                                          max_averkov_N=0,max_heuristic_trials=4);
    if H<>[1,0,2*10^6,0,2*10^6] then error "generator indices were not restored"; end if;
    valid(-1,[1,gens[1],0,gens[2]],H,x);

    # A zero of P with f<=0 blocks every scale, including an included boundary.
    for item in [[x,x^2,cells],[x,x,[[0,1,true,true]]],[-1,0,cells]] do
        result:=StrictlyPositiveCert:-univariate_scale_direction(op(item),x,16,8,time[real]()+10);
        if result<>[FAIL,0] then error "zero-blocked direction was not rejected"; end if;
    end do;
    # Open endpoints do not introduce a spurious zero obstruction.
    result:=StrictlyPositiveCert:-univariate_scale_direction(x,x,[[0,1,false,true]],x,16,8,time[real]()+10);
    if result[1]=FAIL or not StrictlyPositiveCert:-univariate_positive_on_cells(x-result[1]*x,[[0,1,false,true]],x) then
        error "open endpoint direction failed";
    end if;
    # Opposite-sign violations imply no scale; do not keep bisecting.
    result:=StrictlyPositiveCert:-univariate_scale_direction(2*x-1,x,[[-2,-1,true,true],[1,2,true,true]],x,16,8,time[real]()+10);
    if result<>[FAIL,1] then error "infeasible direction was not recognized"; end if;

    result:=StrictlyPositiveCert:-univariate_scale_direction(1/10,x,cells,x,16,8,time[real]()+10);
    if result[1]=FAIL or not type(result[1],rational) or result[1]<=0 or result[1]>=1/10 then
        error "rational halving search failed";
    end if;
    result:=StrictlyPositiveCert:-univariate_scale_direction(1/10,x,cells,x,1,8,time[real]()+10);
    if result<>[FAIL,1] then error "scalar trial limit was ignored"; end if;
    result:=StrictlyPositiveCert:-univariate_scale_direction(1/10,x,cells,x,16,0,time[real]()+10);
    if result<>[FAIL,1] then error "scalar update limit was ignored"; end if;

    # Fixed directions may miss a witness; explicit disable/budgets must still
    # route to Averkov rather than treating a heuristic failure as infeasibility.
    H:=StrictlyPositiveCert:-spCertificates(x^3+1/10,[B,x],x,use_univariate_heuristics=false,max_averkov_N=4);
    valid(x^3+1/10,[B,x],H,x);
    H:=StrictlyPositiveCert:-spCertificates(-1,[-x^2,x^2-1],x,max_heuristic_trials=0,max_averkov_N=16);
    valid(-1,[-x^2,x^2-1],H,x);
    caught:=false;
    try StrictlyPositiveCert:-spCertificates(x^3+1/10,[B,x],x,use_univariate_heuristics=false,max_averkov_N=0);
    catch "Averkov search exhausted max_averkov_N": caught:=true;
    end try;
    if not caught then error "Averkov exponent cap was ignored"; end if;
    caught:=false;
    try StrictlyPositiveCert:-averkov_lemma_7(x,-1,[-x^2,x^2-1/10^6],-x^2,70,max_N=1);
    catch "Averkov search exhausted max_averkov_N": caught:=true;
    end try;
    if not caught then error "Averkov expanded beyond an unsuccessful cap"; end if;

    # The deadline now also limits preprocessing and Averkov, before monogenic.
    caught:=false;
    try StrictlyPositiveCert:-spCertificates(-1,[-x^2,x^2-1/10^6],x,
                       use_univariate_heuristics=false,max_search_time=1/1000);
    catch "certificate search exhausted max_search_time": caught:=true;
    end try;
    if not caught then error "shared public deadline was ignored"; end if;
    for item in [[x+1.5,[B]],[x+1,[B+1.0]]] do
        caught:=false;
        try StrictlyPositiveCert:-spCertificates(op(item),x);
        catch "f and gens must be univariate polynomials with rational coefficients": caught:=true;
        end try;
        if not caught then error "nonrational input was accepted"; end if;
    end do;
$ifdef LOG_TIME
    if StrictlyPositiveCert:-stack_level<>-1 then error "unbalanced heuristic timing stack"; end if;
$endif
    printf("PASS: symbolic division, rational scalar search, exact zero obstructions, index restoration, fallback caps, and shared deadlines\n");
end proc:
check_univariate_heuristics():
