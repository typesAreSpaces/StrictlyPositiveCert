# Compile with LOG_TIME and ENABLE_QUICK_CERT; keep logs outside the checkout.
read "src/StrictlyPositiveCert.mpl":
kernelopts(opaquemodules=false):
if getenv("SPCERT_LOG_DIR")="" then error "SPCERT_LOG_DIR is required"; end if:
currentdir(getenv("SPCERT_LOG_DIR")):

check_logging := proc()
local x, H, item, caught, original_positive;
local balanced := proc(label)
    if StrictlyPositiveCert:-stack_level<>-1 then
        error "unbalanced timing stack after %1",label;
    end if;
end proc;
    H := StrictlyPositiveCert:-spCertificates(x^2+1,[1-x^2],x);
    if H<>[x^2+1,0] then error "nonnegative shortcut changed"; end if;
    balanced("nonnegative shortcut");
    H := StrictlyPositiveCert:-spCertificates(x,[-1],x);
    if expand(x-StrictlyPositiveCert:-dot_product([1,-1],H))<>0 then
        error "negative constant shortcut changed";
    end if;
    balanced("negative constant shortcut");
    H := StrictlyPositiveCert:-spCertificates(x^3+1/10,[1-x^2,x],x,max_averkov_N=0);
    if H<>[1/10,0,x^2] then error "division shortcut changed"; end if;
    balanced("division shortcut");
    H := StrictlyPositiveCert:-spCertificates(-1,[-x^2,x^2-1/10^6],x,max_averkov_N=0);
    if H<>[1,2*10^6,2*10^6] then error "scalar shortcut changed"; end if;
    balanced("scalar shortcut");
    H := StrictlyPositiveCert:-spCertificates(x+3/2,[1,1-x^2,0],x);
    if nops(H)<>4 or H[2]<>0 or H[4]<>0 then error "restoration changed"; end if;
    balanced("generator filtering");
    H := StrictlyPositiveCert:-spCertificates(1-x^2,[2-x^4,1/4-x^2],x);
    if H<>[3/4,0,1] then error "quick certificate shortcut changed"; end if;
    balanced("quick certificate shortcut");
    H := StrictlyPositiveCert:-spCertificates(4/3-x^2,[(1-x^2)^3],x);
    if expand(4/3-x^2-StrictlyPositiveCert:-dot_product([1,(1-x^2)^3],H))<>0 then
        error "full certificate changed";
    end if;
    balanced("full certificate");
    StrictlyPositiveCert:-averkov_lemma_7(x,x+1/2,[1-x^2,x],1-x^2,70);
    balanced("Averkov construction");
    StrictlyPositiveCert:-bound_poly([x+1,1-x],x);
    balanced("odd generator bound");
    if [StrictlyPositiveCert:-bound_poly([x^2+1],x)]<>[] then
        error "bound_poly fall-through changed";
    end if;
    balanced("bound_poly fall-through");
    for item in [[x<=1,x<=2],[x<=1,0<=x],[0<=x,x<=1],[0<=x,1<=x],
                 [x=x],[x=0],[0=x],[x<=1],[0<=x],[x<1]] do
        StrictlyPositiveCert:-bound_info(x,item,0);
        balanced("bound_info branch");
    end do;
    if not StrictlyPositiveCert:-checkPositivityOverSAS([[0<=x,x<=1]],x^2+1,x) then
        error "positivity result changed";
    end if;
    balanced("nested positivity check");

    caught := false;
    try StrictlyPositiveCert:-spCertificates(x,[],x);
    catch "no nonredundant generators remain": caught := true;
    end try;
    if not caught then error "empty-generator error was suppressed"; end if;
    balanced("empty-generator error");
    caught := false;
    try StrictlyPositiveCert:-spCertificates(-1,[1-x^2],x);
    catch "f must be strictly positive": caught := true;
    end try;
    if not caught then error "nested Averkov error was suppressed"; end if;
    balanced("nested Averkov error");
    caught := false;
    try StrictlyPositiveCert:-bound_poly([0,1-x^2],x);
    catch: caught := true;
    end try;
    if not caught then error "bound preprocessing error was suppressed"; end if;
    balanced("bound preprocessing error");
    caught := false;
    try StrictlyPositiveCert:-bound_info(x,[],0);
    catch: caught := true;
    end try;
    if not caught then error "invalid interval error was suppressed"; end if;
    balanced("invalid interval error");
    caught := false;
    original_positive := eval(StrictlyPositiveCert:-univariate_positive_on_cells);
    StrictlyPositiveCert:-univariate_positive_on_cells := proc()
        error "forced positivity failure";
    end proc;
    try StrictlyPositiveCert:-checkPositivityOverSAS([[0<=x,x<=1]],x^2+1,x);
    catch "forced positivity failure": caught := true;
    finally
        StrictlyPositiveCert:-univariate_positive_on_cells := eval(original_positive);
    end try;
    if not caught then error "nested positivity error was suppressed"; end if;
    balanced("positivity error");
    caught := false;
    StrictlyPositiveCert:-univariate_positive_on_cells := proc()
    local spin := proc() while true do end do; end proc;
        timelimit(1/100,spin());
    end proc;
    # The logged parent is above timelimit, so it receives a catchable error.
    # Maple bypasses finally blocks below an outer timelimit boundary.
    try StrictlyPositiveCert:-checkPositivityOverSAS([[0<=x,x<=1]],x^2+1,x);
    catch "time expired": caught := true;
    finally
        StrictlyPositiveCert:-univariate_positive_on_cells := eval(original_positive);
    end try;
    if not caught then error "timeout was suppressed"; end if;
    balanced("catchable timeout");
    printf("PASS: timing cleanup on shortcuts, filtering, normal returns, fall-through, nested errors, and catchable timeouts\n");
end proc:
check_logging():
