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
local averkov_lemma_7 := proc(x, f, gens, B_poly, N_guess)
local cells, gamma := 1, eps, g_i, samples;
local N, N_top := 1, N_bot := 0, N_cur, residual, sample_ok, next_N;
$ifdef LOG_TIME
    INIT_START_LOG_TIME("averkov_lemma_7")
    try
$endif
    cells := univariate_select_cells(univariate_sign_chart(B_poly,x),"nonnegative");
    if univariate_positive_on_cells(f,cells,x) then
        return map(g_i -> 0, gens);
    end if;

    # Bound each positive generator value by 2*gamma on the compact set.
    for g_i in gens do
        while not univariate_positive_on_cells(2*gamma-g_i,cells,x,false) do
            gamma := 2*gamma;
        end do;
    end do;
    eps := findEps(x,gens,f,B_poly);
    DEBUG(__FILE__, __LINE__, lprint(">> Exact Averkov gamma, eps",gamma,eps));

    # A fixed positive prefactor suffices for convergence. Avoid numerical
    # optimization of alpha and logarithmic estimates of the exponent.
    residual := proc(n)
    local term;
        return expand(f-AVERKOV_EXPR(n));
    end proc;
    # A nonpositive value at a point of {B>=0} rejects n without a Sturm
    # sequence; passing samples is necessary, not sufficient.
    samples := averkov_sample_points(cells,B_poly,x);
    sample_ok := proc(n)
    local p := residual(n);
        return andmap(s -> evalb(subs(x=s,p)>0),samples);
    end proc;
    next_N := proc(n)
        # Retain the configured guess as an extra probe when doubling crosses it.
        if n<N_guess and N_guess<2*n then return N_guess; end if;
        return 2*n;
    end proc;
    # Keep N_bot rejected and N_top sample-admissible. Only the final
    # candidate pays for the exact check; if it fails, resume doubling.
    do
        while not sample_ok(N_top) do
            N_bot := N_top;
            N_top := next_N(N_top);
        end do;
$ifdef ENABLE_BINARY_SEARCH_AVKL
        # No monotonicity assumption is needed for correctness; the returned
        # exponent need not be the smallest one.
        while N_top-N_bot>1 do
            N_cur := iquo(N_top+N_bot,2);
            if sample_ok(N_cur) then N_top := N_cur; else N_bot := N_cur; end if;
        end do;
$endif
        if univariate_positive_on_cells(residual(N_top),cells,x) then break; end if;
        N_bot := N_top;
        N_top := next_N(N_top);
    end do;
    N := N_top;
    DEBUG(__FILE__, __LINE__, lprint(">> Exact Averkov exponent",N));
    return map(g_i -> AVERKOV_SOS(g_i,N),gens);
$ifdef LOG_TIME
    finally
        END_LOG_TIME("averkov_lemma_7");
    end try;
$endif
end proc;
