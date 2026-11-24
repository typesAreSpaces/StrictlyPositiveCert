# Rational sample points of {B>=0}: a uniform grid over a rational box
# around each cell, filtered by an exact sign test of B. Used only to reject
# exponents cheaply; acceptance always needs the exact check.
local averkov_sample_points := proc(cells, B, x, per_cell::posint := 64)
local samples := {}, cell, a, b, k, s;
    for cell in cells do
        a := floor(evalf(cell[1])); b := ceil(evalf(cell[2]));
        if a=b then b := a+1; end if;
        for k from 0 to per_cell do
            s := a+k*(b-a)/per_cell;
            if subs(x=s,B)>=0 then samples := samples union {s}; end if;
        end do;
    end do;
    return [op(samples)];
end proc;

# Return SOS multipliers H such that f-dot_product(H,gens)>0 on {B_poly>=0}.
# Every bound and every accepted exponent is checked by exact sign decisions.
local averkov_lemma_7 := proc(x, f, gens, B_poly, N_guess,
                              {gamma_bisection_steps::nonnegint := 4})
local cells, gamma := 1, eps, g_i, samples, values, bounded, gamma_lo, gamma_mid, i;
local N_top := 1, residual, sample_ok, exact_ok, next_N;
$ifdef ENABLE_BINARY_SEARCH_AVKL
local N_bot := 0, N_cur, N_lo, N_hi;
$endif
$ifdef LOG_TIME
    INIT_START_LOG_TIME("averkov_lemma_7")
    try
$endif
    cells := univariate_select_cells(univariate_sign_chart(B_poly,x),"nonnegative");
    if univariate_positive_on_cells(f,cells,x) then
        return map(g_i -> 0, gens);
    end if;

    # Bound each generator value by 2*gamma on the compact set. The exponent
    # grows roughly like gamma/eps, so after doubling, bisect towards the
    # smallest valid bound; every candidate passes the same exact test.
    bounded := c -> andmap(g_i -> univariate_positive_on_cells(2*c-g_i,cells,x,false),gens);
    for g_i in gens do
        while not univariate_positive_on_cells(2*gamma-g_i,cells,x,false) do
            gamma := 2*gamma;
        end do;
    end do;
    if gamma > 1 then
        gamma_lo := gamma/2;
        for i from 1 to gamma_bisection_steps do
            gamma_mid := (gamma_lo+gamma)/2;
            if bounded(gamma_mid) then gamma := gamma_mid; else gamma_lo := gamma_mid; end if;
        end do;
    end if;
    eps := findEps(x,gens,f,B_poly);
    DEBUG(__FILE__, __LINE__, lprint(">> Exact Averkov gamma, eps",gamma,eps));

    # A fixed positive prefactor suffices for convergence. Avoid numerical
    # optimization of alpha and logarithmic estimates of the exponent.
    residual := proc(n)
    local term;
        return expand(f-AVERKOV_EXPR(n));
    end proc;
    # A nonpositive value at a point of {B>=0} rejects n without expanding
    # the residual: [f(s), [g_1(s), ...]] turn each test into rational
    # arithmetic. Passing samples is necessary, not sufficient.
    samples := averkov_sample_points(cells,B_poly,x);
    values := map(s -> [subs(x=s,f), map(g_i -> subs(x=s,g_i),gens)], samples);
    sample_ok := proc(n)
    local v, w;
        for v in values do
            if v[1]-add(AVERKOV_1_EXPR(w,n), w in v[2]) <= 0 then return false; end if;
        end do;
        return true;
    end proc;
    next_N := proc(n)
        # Retain the configured guess as an extra probe when doubling crosses it.
        if n<N_guess and N_guess<2*n then return N_guess; end if;
        return 2*n;
    end proc;
    exact_ok := n -> univariate_positive_on_cells(residual(n),cells,x);
    # Keep N_bot rejected and N_top sample-admissible; resume doubling when
    # N_top fails the exact check. No monotonicity assumption is needed for
    # correctness, and the returned exponent need not be the smallest one.
    do
        while not sample_ok(N_top) do
$ifdef ENABLE_BINARY_SEARCH_AVKL
            N_bot := N_top;
$endif
            N_top := next_N(N_top);
        end do;
$ifdef ENABLE_BINARY_SEARCH_AVKL
        # Bisect (N_bot,N_top] by samples alone. When samples are informative,
        # one exact check accepts the smallest sample-admissible exponent.
        N_lo, N_hi := N_bot, N_top;
        while N_hi-N_lo>1 do
            N_cur := iquo(N_hi+N_lo,2);
            if sample_ok(N_cur) then N_hi := N_cur; else N_lo := N_cur; end if;
        end do;
        if N_hi<N_top then
            if exact_ok(N_hi) then N_top := N_hi; break; end if;
            N_lo := N_hi;
        end if;
$endif
        if exact_ok(N_top) then
$ifdef ENABLE_BINARY_SEARCH_AVKL
            # Samples were not informative below N_top. Bisect between the
            # rejected N_lo and the verified N_top with exact checks; trusting
            # samples here would step the exponent up by one per exact check.
            while N_top-N_lo>1 do
                N_cur := iquo(N_top+N_lo,2);
                if sample_ok(N_cur) and exact_ok(N_cur) then N_top := N_cur; else N_lo := N_cur; end if;
            end do;
$endif
            break;
        end if;
$ifdef ENABLE_BINARY_SEARCH_AVKL
        N_bot := N_top;
$endif
        N_top := next_N(N_top);
    end do;
    DEBUG(__FILE__, __LINE__, lprint(">> Exact Averkov exponent",N_top));
    return map(g_i -> AVERKOV_SOS(g_i,N_top),gens);
$ifdef LOG_TIME
    finally
        END_LOG_TIME("averkov_lemma_7");
    end try;
$endif
end proc;
