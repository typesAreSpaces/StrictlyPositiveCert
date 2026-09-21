# We assume:
# 1. SemiAlgebraic(g) is bounded
# Return: a polynomial l
# such that poly - l has a lowerbound
# over \mathbb{R}
# and an epsilon (eps_LS) for Last_step
local lower_bound_poly := proc(x, f, g)
$ifdef LOG_TIME
    INIT_START_LOG_TIME("lower_bound_poly",0)
$endif
local d_f, c_f;
local d_g, d_diff;
local A, disc, roots_disc;
local h, _point, _point_candidates;
local c;
# The following is used in the
# minimization problem to find
# c in order to avoid the boundary
# points
local eps := 1/100;
# This variable is passed to Last_step
local eps_LS := -1, curr_eps_LS := -1;

    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> f", f));
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> g", g));
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
        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> f is bounded: "));
$ifdef LOG_TIME
        END_LOG_TIME("lower_bound_poly",0)
$endif
        return 0, eps_LS;
    end if;

    d_g := degree(expand(g), x); # quick_degree
    d_diff := d_f - d_g;
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> d_diff", d_diff));
    if 0 <= d_diff then
        if type(d_diff, even) then
            d_diff := d_diff + 2;
        else
            d_diff := d_diff + 1;
        end if;
    else
        c := find_constant_lower_bound_poly(f, 1, g, x, eps);
        return c, eps_LS;
    end if;

    disc := diff(f,x)*g - f*diff(g, x);
    A := x - d_diff*f*g/disc;
    roots_disc := select(_root -> evalf(subs(_root, g)) > 0,
                         Isolate(disc, maxprec=1000, digits=30));
    _point_candidates := map(
        _root -> convert(subs({x=round(op(_root)[2]*1000)/1000}, A), rational),
        roots_disc);

# Loop to choose optimal _point
    for _point in _point_candidates do
        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> _point", _point));

# TODO Compute h using 'more diverse' _points
        h := (x - _point)^d_diff;
        DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> h", h));

        c := find_constant_lower_bound_poly(f, h, g, x, eps);

# We want is maximize eps_LS
$ifdef WEIFENG_OPTIMIZATION
        curr_eps_LS := evalf(findEps(x, [g], f));
$else
        curr_eps_LS := evalf(findEps(x, [g], f - c*h*g));
$endif
        curr_eps_LS := 1/2*convert(curr_eps_LS, rational);
        if eps_LS < curr_eps_LS then
            eps_LS := curr_eps_LS;
            DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">>> evalf(_point)", evalf(_point)));
            DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">>> _point", _point));
            DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">>> current eps_LS", evalf(eps_LS)));
        end if;
    end do;
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">>> Final eps_LS", evalf(eps_LS)));
    return c*h, eps_LS;
end proc;

local find_constant_lower_bound_poly := proc(f, h, g, x, eps)
local c, G := h*g;
local opt_roots := Isolate(diff(f,x)*G - f*diff(G, x), maxprec=1000, digits=30);
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> opt_roots", opt_roots));

# We just need a lowerbound, not the
# tightest lowerbound [to discuss later]
# TODO Figure out `optimal' constant (i.e., 9/10, 99/100, ...)
# to avoid eps_LS become a negative number
#c := 9/10*min(map(x_arg -> subs(x_arg, f/G), select(_root-> evalf(subs(_root, g)) > 0, opt_roots)));
    c := 999/1000*min(
        map(x_arg ->
            if (evalf(subs(x_arg, G)) < eps) then
                # If we are minimizing over
                # an isolated point of S(g), any value of
                # C satisfy the minimization condition
                1
            else
                subs(x_arg, f/G)
            end if,
            select(_root-> evalf(subs(_root, g)) >= 0, opt_roots))
                     );

    c := convert(evalf(c), rational);
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> c as rational", c));

    return c;
end proc;
