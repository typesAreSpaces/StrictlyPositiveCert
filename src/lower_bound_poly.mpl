# We assume:
# 1. SemiAlgebraic(g) is bounded (so g has even degree and lcoeff(g) < 0)
# 2. f > 0 over SemiAlgebraic(g)
# Return: a sum of squares s = c*q such that
# f - s*g has a lowerbound over \mathbb{R}
# and is positive over SemiAlgebraic(g),
# and an epsilon (eps_LS) for Last_step.
#
# Whenever possible, c is chosen so that f - s*g is
# strictly positive over \mathbb{R}; then Last_step
# returns 0 immediately and eps_LS = -1 is returned.
local lower_bound_poly := proc(x, f, g)
$ifdef LOG_TIME
    INIT_START_LOG_TIME("lower_bound_poly",0)
$endif
local d_f, c_f;
local d_g, d_diff;
local k, k_first, k_bounded;
local L, shapes, q, c, c_low, c_high, margin;
local best_c, best_q, best_margin;
# This variable is passed to Last_step
local eps_LS := -1, curr_eps_LS;

    DEBUG(__FILE__, __LINE__, lprint(">> f", f));
    DEBUG(__FILE__, __LINE__, lprint(">> g", g));
$ifdef LOG_TIME
    START_LOG_TIME("lower_bound_poly::expand(f)",1);
$endif
    d_f := degree(expand(f), x); # quick_degree
    c_f := coeff(f, x^d_f);
$ifdef LOG_TIME
    END_LOG_TIME("lower_bound_poly::expand(f)",1);
$endif

# If f has a lowerbound over \mathbb{R}
# then make no changes to f
    if type(d_f, even) and evalb(evalf(c_f) > 0) then
        DEBUG(__FILE__, __LINE__, lprint(">> f is bounded: "));
$ifdef LOG_TIME
        END_LOG_TIME("lower_bound_poly",0)
$endif
        return 0, eps_LS;
    end if;

    d_g := degree(expand(g), x); # quick_degree
    d_diff := d_f - d_g;
    DEBUG(__FILE__, __LINE__, lprint(">> d_diff", d_diff));
# deg(q) = 2*k. With k >= k_bounded, f - c*q*g is bounded
# below for every c > 0. With k_first <= k < k_bounded,
# the leading coefficient depends on c; ratio_bounds
# accounts for it through the limit at infinity.
    k_bounded := `if`(d_diff < 0, 0, floor(d_diff/2) + 1);
    k_first := max(0, ceil(d_diff/2));
# The hull is only needed for shapes of positive degree
    L := FAIL;

# Look for the smallest degree where f - c*q*g is
# strictly positive over \mathbb{R}
    for k from k_first to k_bounded + 1 do
        best_q := 0;
        if k > 0 and L = FAIL then
            L := hull_affine_map(g, x);
            DEBUG(__FILE__, __LINE__, lprint(">> L", L));
        end if;
        for q in lbp_shapes(k, L, x) do
            c_low, c_high := ratio_bounds(f, q, g, x);
            DEBUG(__FILE__, __LINE__, lprint(">> k, c_low, c_high", k, c_low, c_high));
            if not evalb(evalf(c_low) < evalf(c_high)) then
                next;
            end if;
            if c_high = infinity then
                c := `if`(c_low = 0, 1, 2*c_low);
                margin := infinity;
            elif c_low = 0 then
                c := c_high/2;
                margin := infinity;
            else
                c := sqrt(c_low*c_high);
                margin := c_high/c_low;
            end if;
            c := convert(evalf(c), rational, 6);
# The bounds are numerical; the certificate is checked exactly
            if isSOS(expand(f - c*q*g), x) and
               (best_q = 0 or evalf(best_margin < margin)) then
                best_c, best_q, best_margin := c, q, margin;
            end if;
        end do;
        if best_q <> 0 then
            DEBUG(__FILE__, __LINE__, lprint(">> f - s*g is SOS with s", best_c*best_q));
$ifdef LOG_TIME
            END_LOG_TIME("lower_bound_poly",0)
$endif
            return best_c*best_q, -1;
        end if;
    end do;

# Otherwise keep f - c*q*g positive over SemiAlgebraic(g)
# with the guaranteed degree, maximizing eps_LS.
# Candidate, constant and eps_LS are kept together.
    if k_bounded > 0 and L = FAIL then
        L := hull_affine_map(g, x);
    end if;
    shapes := lbp_shapes(k_bounded, L, x);
    if nops(shapes) = 1 then
# Nothing to choose: Last_step computes the same eps itself
$ifdef LOG_TIME
        END_LOG_TIME("lower_bound_poly",0)
$endif
        return find_constant_lower_bound_poly(f, shapes[1], g, x)*shapes[1], eps_LS;
    end if;
    best_q := 0;
    for q in shapes do
        c := find_constant_lower_bound_poly(f, q, g, x);
$ifdef WEIFENG_OPTIMIZATION
        curr_eps_LS := evalf(findEps(x, [g], f));
$else
        curr_eps_LS := evalf(findEps(x, [g], f - c*q*g));
$endif
        curr_eps_LS := 1/2*convert(curr_eps_LS, rational);
        if best_q = 0 or eps_LS < curr_eps_LS then
            best_c, best_q, eps_LS := c, q, curr_eps_LS;
            DEBUG(__FILE__, __LINE__, lprint(">>> current eps_LS", evalf(eps_LS)));
        end if;
    end do;
    DEBUG(__FILE__, __LINE__, lprint(">>> Final eps_LS", evalf(eps_LS)));
$ifdef LOG_TIME
    END_LOG_TIME("lower_bound_poly",0)
$endif
    return best_c*best_q, eps_LS;
end proc;

# Returns L(x) = (x - m)/r mapping (a rational outer
# approximation of) the convex hull of SemiAlgebraic(g)
# onto [-1, 1]. Since g < 0 near infinity, the hull is
# spanned by the smallest and largest real roots of g.
local hull_affine_map := proc(g, x)
local roots_g := map(rhs, Isolate(g, x, digits=30));
local left, right, r;
    if roots_g = [] then
        return x;
    end if;
    left := floor(1000*min(roots_g))/1000;
    right := ceil(1000*max(roots_g))/1000;
    r := (right - left)/2;
    if r = 0 then
        r := 1;
    end if;
    return (x - (left + right)/2)/r;
end proc;

# Candidate squares of degree 2*k: the centered monomial
# L^(2k) and the Chebyshev square T_k(L)^2. Both are at most
# 1 over the hull; the latter grows faster right outside it.
local lbp_shapes := proc(k, L, x)
    if k = 0 then
        return [1];
    elif k = 1 then
        return [expand(L^2)];
    end if;
    return [expand(L^(2*k)), expand(orthopoly[T](k, L)^2)];
end proc;

# f - c*q*g >= 0 over \mathbb{R} iff c_low <= c <= c_high, where
#   c_low  = sup of f/(q*g) over {g < 0},
#   c_high = inf of f/(q*g) over {g > 0}.
# Near zeros of g the ratio goes to -infinity outside and to
# +infinity inside SemiAlgebraic(g), so only critical points,
# zeros of q, and the limit at infinity matter.
local ratio_bounds := proc(f, q, g, x)
local G := expand(q*g);
local d_f := degree(expand(f), x), d_G := degree(G, x);
local c_low := 0, c_high := infinity;
local num, _root, p, g_p, G_p;
    num := expand(diff(f, x)*G - f*diff(G, x));
    if num = 0 then
        p := normal(f/G);
        return p, p;
    end if;
# SemiAlgebraic(g) is bounded, so g < 0 near infinity
    if d_G < d_f then
        c_low := infinity;
    elif d_G = d_f then
        c_low := lcoeff(expand(f), x)/lcoeff(G, x);
    end if;
# A zero of q where g < 0 and f <= 0 is a pole towards +infinity
    if degree(q, x) > 0 then
        for _root in Isolate(q, x, digits=30) do
            p := rhs(_root);
            if evalf(subs(x=p, g)) < 0 and evalf(subs(x=p, f)) <= 0 then
                c_low := infinity;
            end if;
        end do;
    end if;
    for _root in Isolate(num, x, maxprec=1000, digits=30) do
        p := rhs(_root);
        g_p := evalf[30](subs(x=p, g));
        G_p := evalf[30](subs(x=p, G));
        if G_p = 0 then
            next;
        end if;
        if g_p > 0 then
            c_high := min(c_high, evalf[30](subs(x=p, f))/G_p);
        elif g_p < 0 then
            c_low := max(c_low, evalf[30](subs(x=p, f))/G_p);
        end if;
    end do;
    return c_low, c_high;
end proc;

# Largest (up to 999/1000) constant c such that
# f - c*h*g > 0 over SemiAlgebraic(g)
local find_constant_lower_bound_poly := proc(f, h, g, x)
local c_low, c_high;
    c_low, c_high := ratio_bounds(f, h, g, x);
    DEBUG(__FILE__, __LINE__, lprint(">> c_high", c_high));
# No critical point inside SemiAlgebraic(g) (e.g., finitely
# many points): h*g vanishes there and any c > 0 works
    if c_high = infinity then
        return 1;
    end if;
    return convert(evalf(999/1000*c_high), rational);
end proc;
