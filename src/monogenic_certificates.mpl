# Exact monogenic certificates, ported from monogenicCert.wl.
# Return [sigma0, sigma1] with f = sigma0 + sigma1*g.  The default output
# contains SOS polynomials, matching spCertificates.  With output_sos=true,
# each sigma is a list of [weight, polynomial] pairs representing weighted squares.
# The optional decomposition uses Schweighofer's quadratic reduction (univsos1).
local monogenic_positive_on_R := proc(poly, x)
local p := expand(poly), d;
    if p = 0 then return false; end if;
    d := degree(p, x);
    return type(d, even) and lcoeff(p, x) > 0
    and (d = 0 or sturm(p, x, -infinity, infinity) = 0);
end proc;

local monogenic_status_lambda := proc(f, g, sg, lambda, x, regions := FAIL)
local h := expand(f-lambda*sg), too_big, too_small, cells := regions;
    if monogenic_positive_on_R(h, x) then return "ok"; end if;
    if cells=FAIL then cells := monogenic_sign_regions(g,x); end if;
    too_big := not univariate_positive_on_cells(h,cells[1],x);
    too_small := not univariate_positive_on_cells(h,cells[3],x);
    if too_big and too_small then return "infeasible";
    elif too_big then return "big";
    else return "small";
    end if;
end proc;

# One resumable search trial. Polynomial construction is included in its time limit.
local monogenic_search_trial := proc(f, g, x, k, M, base0, sg0, lambda, regions := FAIL)
local base := base0, sg := sg0;
    if base = FAIL then
        base := expand(((M-g)/M)^k);
        sg := expand(base^2*g);
    end if;
    return [base, sg, monogenic_status_lambda(f, g, sg, lambda, x, regions)];
end proc;

# The sign chart of g is shared by normalization and every lambda trial.
local monogenic_sign_regions := proc(g, x)
local chart := univariate_sign_chart(g,x);
    return [univariate_select_cells(chart,"positive"),
            univariate_select_cells(chart,"nonnegative"),
            univariate_select_cells(chart,"nonpositive")];
end proc;

# Normalize once, under the overall search time limit. Keep the original
# dyadic choices of m and M, including equality and the lower cap for M.
local monogenic_search_bounds := proc(f, g, x)
local m := 1, M := 1, regions := monogenic_sign_regions(g,x);
    if regions[2] <> [] then
        if not univariate_positive_on_cells(f,regions[2],x) then
            error "f must be strictly positive on the nonnegative set of g";
        end if;
        while univariate_positive_on_cells(f-2*m,regions[2],x,false) do m:=2*m; end do;
        while not univariate_positive_on_cells(f-m,regions[2],x,false) do m:=m/2; end do;
        # Since M>0, {g>=0,g>M} equals {g>M} on the whole real line.
        while not univariate_nonnegative_on_R(M-g,x) do M:=2*M; end do;
        while M>1/2^20 and univariate_nonnegative_on_R(M/2-g,x) do M:=M/2; end do;
    end if;
    return [m,M,regions];
end proc;

# Bisect a rational isolating interval for a square-free polynomial.
local monogenic_refine_interval := proc(q, x, interval, width)
local a := interval[1], b := interval[2], midpoint, qa, qm;
    while b-a > width do
        midpoint := (a+b)/2;
        qm := subs(x=midpoint, q);
        if qm = 0 then return [midpoint, midpoint]; end if;
        qa := subs(x=a, q);
        if qa = 0 then return [a, a]; end if;
        if subs(x=b, q) = 0 then return [b, b]; end if;
        if qa*qm < 0 then b := midpoint; else a := midpoint; end if;
    end do;
    return [a, b];
end proc;

local monogenic_sos_step := proc(f, x, max_sos_iterations)
local n := degree(f, x), a, b, e, fp, q, intervals, width, candidates;
local iteration, interval, t, ft, fpt, p, r;
    if n = 0 then return [[[f, 1]], []]; end if;
    if n = 2 then
        a := coeff(f, x, 2); b := coeff(f, x, 1); e := coeff(f, x, 0);
        return [[[a, x+b/(2*a)], [e-b^2/(4*a), 1]], []];
    end if;
    fp := diff(f, x);
    q := quo(fp, gcd(fp, diff(fp, x)), x);
    intervals := realroot(q);
    for iteration from 1 to max_sos_iterations do
        width := 1/2^iteration;
        intervals := map(interval -> monogenic_refine_interval(q, x, interval, width), intervals);
        candidates := {};
        for interval in intervals do
            candidates := candidates union
            {interval[1]-width, interval[1], interval[2], interval[2]+width};
        end do;
        for t in sort(convert(candidates, list)) do
            ft := subs(x=t, f); fpt := subs(x=t, fp);
            p := quo(expand(f-ft-fpt*(x-t)), (x-t)^2, x);
            r := expand(p-fpt^2/(4*ft));
            if monogenic_positive_on_R(r, x) then
                return [[[ft, 1+fpt/(2*ft)*(x-t)]], [t, r]];
            end if;
        end do;
    end do;
    return FAIL;
end proc;

local monogenic_exact_sos := proc(poly, x, max_sos_iterations)
local f := expand(poly), squares := [], pref := 1, step, term;
    if f = 0 then return []; end if;
    if not monogenic_positive_on_R(f, x) then
        error "SOS decomposition requires a strictly positive polynomial";
    end if;
    while true do
        step := monogenic_sos_step(f, x, max_sos_iterations);
        if step = FAIL then
            error "exact SOS decomposition exhausted max_sos_iterations=%1", max_sos_iterations;
        end if;
        for term in step[1] do
            squares := [op(squares), [term[1], expand(pref*term[2])]];
        end do;
        if step[2] = [] then return squares; end if;
        pref := pref*(x-step[2][1]);
        f := step[2][2];
    end do;
end proc;

monogenicCertificates := proc(f0, g0, x::name,
                              {max_k::nonnegint := 100, max_bisection::nonnegint := 200,
                               query_time_limit::positive := 1/10, max_query_retries::nonnegint := 4,
                               max_search_time::positive := 300,
                               max_sos_iterations::nonnegint := 64, output_sos::truefalse := false})
local f := expand(f0), g := expand(g0), dg, df, m, M, bounds, regions;
local k, next_k := 0, round := 0, active := 0, lambda, lo, hi, iteration, retries;
local states := table(), bases := table(), products := table();
local trial, status, base, h, started := time[real](), remaining, allowance;
    if not type(f, polynom(rational, x)) or not type(g, polynom(rational, x)) then
        error "f and g must be univariate polynomials with rational coefficients";
    end if;
    if g = 0 then error "g must have a compact nonnegative set"; end if;
    dg := degree(g, x); df := degree(f, x);
    if not type(dg, even) or lcoeff(g, x) >= 0 then
        error "g must have even degree and negative leading coefficient";
    end if;
    bounds := timelimit(max_search_time, monogenic_search_bounds(f, g, x));
    m, M, regions := op(bounds);

    # Admit higher degrees gradually; advance every active lambda search once.
    # Square-spaced admission avoids building many high-degree states while a
    # low-degree search is still making progress. If none remain, admit immediately.
    # State = [lambda, lower bracket, upper bracket, completed trials, retries].
    while active > 0 or next_k <= max_k do
        if next_k <= max_k and (round >= next_k^2 or active = 0) then
            k := next_k;
            next_k := next_k+1;
            states[k] := FAIL;
            if not (dg*(2*k+1) < df and
                    not (type(df, even) and lcoeff(f, x) > 0)) then
                states[k] := [m*(2*k+1)/(2*M), 0, 0, 0, 0];
                bases[k] := FAIL;
                products[k] := FAIL;
                active := active+1;
            end if;
        end if;
        for k from 0 to next_k-1 do
            if states[k] = FAIL then next; end if;
            lambda, lo, hi, iteration, retries := op(states[k]);
            remaining := max_search_time-(time[real]()-started);
            if remaining <= 0 then
                error "monogenic certificate search exhausted max_search_time=%1", max_search_time;
            end if;
            allowance := min(query_time_limit*2^retries, remaining);
            try
                trial := timelimit(allowance,
                                   monogenic_search_trial(f, g, x, k, M, bases[k], products[k], lambda, regions));
            catch "time expired":
                DEBUG(__FILE__, __LINE__, lprint(">> Deferred monogenic trial k, iteration, retry", k, iteration, retries));
                if retries >= max_query_retries then
                    states[k] := FAIL;
                    bases[k] := FAIL; products[k] := FAIL;
                    active := active-1;
                else
                    # Keep the same lambda and brackets; timeouts prove nothing.
                    states[k] := [lambda, lo, hi, iteration, retries+1];
                end if;
                next;
            end try;
            bases[k], products[k], status := op(trial);
            if status = "ok" then
                base := bases[k];
                h := expand(f-lambda*products[k]);
                DEBUG(__FILE__, __LINE__, lprint(">> monogenic certificate k, lambda", k, lambda));
                if output_sos then
                    return [monogenic_exact_sos(h, x, max_sos_iterations), [[lambda, base]]];
                end if;
                return [h, expand(lambda*base^2)];
            end if;
            if status = "infeasible" or iteration >= max_bisection then
                states[k] := FAIL;
                bases[k] := FAIL; products[k] := FAIL;
                active := active-1;
                next;
            end if;
            if status = "small" then lo := lambda; else hi := lambda; end if;
            if lo = 0 then lambda := hi/2;
            elif hi = 0 then lambda := 2*lo;
            else lambda := (lo+hi)/2;
            end if;
            states[k] := [lambda, lo, hi, iteration+1, 0];
        end do;
        round := round+1;
    end do;
    error "monogenic certificate search exhausted max_k=%1, max_bisection=%2, and max_query_retries=%3", max_k, max_bisection, max_query_retries;
end proc;
