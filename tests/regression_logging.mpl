# Compile with LOG_TIME and ENABLE_QUICK_CERT; keep logs outside the checkout.
read "src/StrictlyPositiveCert.mpl":
kernelopts(opaquemodules=false):
if getenv("SPCERT_LOG_DIR")="" then error "SPCERT_LOG_DIR is required"; end if:
currentdir(getenv("SPCERT_LOG_DIR")):

check_logging := proc()
local x, H, item, caught;
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
    H := StrictlyPositiveCert:-spCertificates(x+3/2,[1,1-x^2,0],x);
    if nops(H)<>4 or H[2]<>0 or H[4]<>0 then error "restoration changed"; end if;
    balanced("recursive generator filtering");
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
    for item in [[x^2+1], [x^3,x], [-x^3,-x], []] do
        caught := false;
        try StrictlyPositiveCert:-bound_poly(item,x);
        catch "no compactness witness found": caught := true;
        end try;
        if not caught then error "missing compactness witness was accepted: %1",item; end if;
        balanced("bound_poly compactness error");
    end do;
    caught := false;
    try StrictlyPositiveCert:-spCertificates(x,[x^2+1],x);
    catch "bound_poly could not witness compactness of gens": caught := true;
    end try;
    if not caught then error "public compactness error was suppressed"; end if;
    balanced("public compactness error");
    for item in [[x<=1,x<=2],[x<=1,0<=x],[0<=x,x<=1],[0<=x,1<=x],
                 [x=x],[x=0],[0=x],[x<=1],[0<=x],[x<1]] do
        StrictlyPositiveCert:-bound_info(x,item,0);
        balanced("bound_info branch");
    end do;
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
    printf("PASS: timing cleanup on shortcuts, recursion, normal returns, compactness errors, and nested errors\n");
end proc:
check_logging():
