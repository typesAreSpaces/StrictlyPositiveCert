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

local monogenic_status_lambda := proc(f, g, sg, lambda, x)
local h := expand(f-lambda*sg), too_big, too_small;
    if monogenic_positive_on_R(h, x) then return "ok"; end if;
    too_big := evalb(SemiAlgebraic([g > 0, h <= 0], [x]) <> []);
    too_small := evalb(SemiAlgebraic([g <= 0, h <= 0], [x]) <> []);
    if too_big and too_small then return "infeasible";
    elif too_big then return "big";
    else return "small";
    end if;
end proc;

local monogenic_find_lambda := proc(f, g, sg, lambda0, x, max_bisection)
local lambda := lambda0, lo := 0, hi := 0, status, iteration;
    for iteration from 0 to max_bisection do
        status := monogenic_status_lambda(f, g, sg, lambda, x);
        if status = "ok" then return lambda;
        elif status = "infeasible" then return FAIL;
        elif status = "small" then lo := lambda;
        else hi := lambda;
        end if;
        if lo = 0 then lambda := hi/2;
        elif hi = 0 then lambda := 2*lo;
        else lambda := (lo+hi)/2;
        end if;
    end do;
    return FAIL;
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
     max_sos_iterations::nonnegint := 64, output_sos::truefalse := false})
local f := expand(f0), g := expand(g0), dg, df, m := 1, M := 1;
local k, base, sg, lambda, h, squares;
    if not type(f, polynom(rational, x)) or not type(g, polynom(rational, x)) then
        error "f and g must be univariate polynomials with rational coefficients";
    end if;
    if g = 0 then error "g must have a compact nonnegative set"; end if;
    dg := degree(g, x); df := degree(f, x);
    if not type(dg, even) or lcoeff(g, x) >= 0 then
        error "g must have even degree and negative leading coefficient";
    end if;
    if SemiAlgebraic([g >= 0], [x]) <> [] then
        if SemiAlgebraic([g >= 0, f <= 0], [x]) <> [] then
            error "f must be strictly positive on the nonnegative set of g";
        end if;
        while SemiAlgebraic([g >= 0, f < 2*m], [x]) = [] do m := 2*m; end do;
        while SemiAlgebraic([g >= 0, f < m], [x]) <> [] do m := m/2; end do;
        while SemiAlgebraic([g >= 0, g > M], [x]) <> [] do M := 2*M; end do;
        while M > 1/2^20 and SemiAlgebraic([g >= 0, g > M/2], [x]) = [] do
            M := M/2;
        end do;
    end if;
    for k from 0 to max_k do
        base := expand(((M-g)/M)^k);
        sg := expand(base^2*g);
        if degree(sg, x) < df and not (type(df, even) and lcoeff(f, x) > 0) then
            next;
        end if;
        lambda := monogenic_find_lambda(f, g, sg, m*(2*k+1)/(2*M), x, max_bisection);
        if lambda <> FAIL then
            h := expand(f-lambda*sg);
            DEBUG(__FILE__, __LINE__, lprint(">> monogenic certificate k, lambda", k, lambda));
            if output_sos then
                squares := monogenic_exact_sos(h, x, max_sos_iterations);
                return [squares, [[lambda, base]]];
            end if;
            return [h, expand(lambda*base^2)];
        end if;
    end do;
    error "monogenic certificate search exhausted max_k=%1 and max_bisection=%2", max_k, max_bisection;
end proc;
