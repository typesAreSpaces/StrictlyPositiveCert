# Check scheduler fairness, retry state, limits, and propagation of solver errors.
read "src/StrictlyPositiveCert.mpl":
kernelopts(opaquemodules=false):

check_scheduler := proc()
local x, f := 4/3-x^2, g := (1-x^2)^3, cert, original_trial;
local visits := [], caught;
    original_trial := eval(StrictlyPositiveCert:-monogenic_search_trial);
    try
        # Simulate a slow k=0 query while allowing the real k=1 search to finish.
        StrictlyPositiveCert:-monogenic_search_trial := proc(f, g, x, k, M, base, sg, lambda)
            visits := [op(visits), [k, lambda]];
            if k = 0 then error "time expired"; end if;
            return original_trial(f, g, x, k, M, base, sg, lambda);
        end proc;
        cert := StrictlyPositiveCert:-spCertificates(f, [g], x,
            max_k=1, max_query_retries=1, max_search_time=10);
        if expand(f-cert[1]-cert[2]*g) <> 0
           or expand(cert[2]-3/2*(1-g)^2) <> 0 then
            error "deferred query blocked the later exponent";
        end if;
        if visits <> [[0,1/8], [0,1/8], [1,3/8], [1,3/4], [1,3/2]] then
            error "retry changed a bracket, restarted a search, or ignored its limit: %1", visits;
        end if;
        visits := [];
        caught := false;
        try
            StrictlyPositiveCert:-spCertificates(f,[g],x,
                max_k=0,max_query_retries=0);
        catch: caught := evalb(StringTools:-Search("search exhausted",sprintf("%a",[lastexception])) > 0);
        end try;
        if not caught or visits <> [[0,1/8]] then
            error "public retry budget was ignored";
        end if;

        # A timed-out candidate can succeed on retry without consuming bisections.
        visits := [];
        StrictlyPositiveCert:-monogenic_search_trial := proc(f, g, x, k, M, base, sg, lambda)
            visits := [op(visits), [k, lambda]];
            if nops(visits) = 1 then error "time expired"; end if;
            return original_trial(f, g, x, k, M, base, sg, lambda);
        end proc;
        cert := StrictlyPositiveCert:-monogenicCertificates(x+3/2,1-x^2,x,
            max_k=0,max_bisection=0,max_query_retries=1);
        if visits <> [[0,1/4],[0,1/4]]
           or expand(x+3/2-cert[1]-cert[2]*(1-x^2)) <> 0 then
            error "transient timeout lost its candidate or consumed its trial budget";
        end if;

        # Errors other than timeouts must reach the caller unchanged.
        StrictlyPositiveCert:-monogenic_search_trial := proc()
            error "scheduler test solver failure";
        end proc;
        caught := false;
        try StrictlyPositiveCert:-monogenicCertificates(f,g,x);
        catch "scheduler test solver failure": caught := true;
        end try;
        if not caught then error "solver error was hidden as a timeout"; end if;

        # Exercise a real cooperative query timeout, not just the mock above.
        StrictlyPositiveCert:-monogenic_search_trial := proc()
        local j := 0;
            while true do j := j+1; end do;
        end proc;
        caught := false;
        try
            StrictlyPositiveCert:-monogenicCertificates(f,g,x,
                max_k=0,query_time_limit=1/100,max_query_retries=0);
        catch: caught := evalb(StringTools:-Search("search exhausted",sprintf("%a",[lastexception])) > 0);
        end try;
        if not caught then error "query timeout did not retire its exhausted state"; end if;

        # A deadline must stop retries before all exponent/trial budgets run out.
        caught := false;
        try
            StrictlyPositiveCert:-spCertificates(f,[g],x,
                max_search_time=1/10,query_time_limit=1/100,max_query_retries=20);
        catch: caught := evalb(StringTools:-Search("max_search_time",sprintf("%a",[lastexception])) > 0);
        end try;
        if not caught then error "public search deadline was ignored"; end if;
    finally
        StrictlyPositiveCert:-monogenic_search_trial := original_trial;
    end try;
    cert := timelimit(10,StrictlyPositiveCert:-spCertificates(f,[g],x));
    if expand(f-cert[1]-cert[2]*g) <> 0 or degree(cert[1],x) <> 18
       or SemiAlgebraic([cert[1]<=0],[x]) <> [] then
        error "default scheduler did not resolve the Stengle regression";
    end if;
    printf("PASS: interleaved search, retry state, query limits, deadline, solver errors, and Stengle regression\n");
end proc:
check_scheduler():
