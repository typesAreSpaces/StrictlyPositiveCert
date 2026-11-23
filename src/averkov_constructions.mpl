# Return SOS multipliers H such that f-dot_product(H,gens)>0 on {B_poly>=0}.
# Every bound and every accepted exponent is checked by exact sign decisions.
# Time is bounded by the caller's timelimit; max_N bounds the exponent search.
local averkov_lemma_7 := proc(x, f, gens, B_poly, N_guess, {max_N::nonnegint := 256})
local cells, gamma := 1, eps, g_i;
local N, N_top := 1, next_N, admissible;
$ifdef ENABLE_BINARY_SEARCH_AVKL
local N_bot := 0, N_cur;
$endif
$ifdef LOG_TIME
    INIT_START_LOG_TIME("averkov_lemma_7",0)
    try
$endif
    cells := univariate_select_cells(univariate_sign_chart(B_poly,x),"nonnegative");
    if univariate_positive_on_cells(f,cells,x) then
        return map(g_i -> 0, gens);
    end if;

    if max_N=0 then error "Averkov search exhausted max_averkov_N=%1",max_N; end if;

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
    admissible := proc(n)
    local term;
        return univariate_positive_on_cells(expand(f-AVERKOV_EXPR(n)),cells,x);
    end proc;
    while not admissible(N_top) do
$ifdef ENABLE_BINARY_SEARCH_AVKL
        N_bot := N_top;
$endif
        if N_top>=max_N then
            error "Averkov search exhausted max_averkov_N=%1",max_N;
        end if;
        next_N := min(2*N_top,max_N);
        # Retain the configured guess as an extra probe within the cap.
        if N_top<N_guess and N_guess<next_N then next_N:=N_guess; end if;
        N_top := next_N;
    end do;

$ifdef ENABLE_BINARY_SEARCH_AVKL
    # Keep a verified upper endpoint. No monotonicity assumption is needed
    # for correctness; the returned exponent need not be the smallest one.
    while N_top-N_bot>1 do
        N_cur := iquo(N_top+N_bot,2);
        if admissible(N_cur) then
            N_top := N_cur;
        else
            N_bot := N_cur;
        end if;
    end do;
$endif
    N := N_top;
    DEBUG(__FILE__, __LINE__, lprint(">> Exact Averkov exponent",N));
    return map(g_i -> AVERKOV_SOS(g_i,N),gens);
$ifdef LOG_TIME
    finally
        END_LOG_TIME("averkov_lemma_7",0);
    end try;
$endif
end proc;
