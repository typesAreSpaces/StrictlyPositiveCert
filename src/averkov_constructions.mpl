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

# Chebyshev polynomial T_n(u), expanded, for a rational or polynomial u.
# The ladder keeps (T_k, T_{k+1}) and doubles k bit by bit, using
#   T_{2k} = 2*T_k^2-1,  T_{2k+1} = 2*T_k*T_{k+1}-u,  T_{2k+2} = 2*T_{k+1}^2-1,
# so it needs O(log n) products instead of the n of the three-term recurrence.
local averkov_chebyshev := proc(n::nonnegint, u0)
local u := expand(u0), a := 1, b, bit;
    if n=0 then return 1; end if;
    b := u;
    for bit in ListTools:-Reverse(convert(n,base,2)) do
        if bit=1 then
            a, b := expand(2*a*b-u), expand(2*b^2-1);
        else
            a, b := expand(2*a^2-1), expand(2*a*b-u);
        end if;
    end do;
    return a;
end proc;

# A large rational c with p-c>0 on cells: double while valid (or halve until
# valid), then bisect. Cells are compact and p is strictly positive on them,
# so both loops stop; the caps are safety nets. Empty cells impose nothing.
local averkov_positive_lower_bound := proc(p, cells, x,
                                          {bisection_steps::nonnegint := 3, max_steps::posint := 200})
local ok, lo, hi, mid, i;
    if cells=[] then return 1; end if;
    ok := c -> univariate_positive_on_cells(p-c,cells,x,true);
    if ok(1) then
        lo, hi := 1, 2;
        for i from 1 to max_steps while ok(hi) do lo, hi := hi, 2*hi; end do;
    else
        hi := 1;
        for i from 1 to max_steps while not ok(hi/2) do hi := hi/2; end do;
        lo := hi/2;
        if not ok(lo) then error "no positive lower bound found on the separated set"; end if;
    end if;
    for i from 1 to bisection_steps do
        mid := (lo+hi)/2;
        if ok(mid) then lo := mid; else hi := mid; end if;
    end do;
    return lo;
end proc;

# A rational F>=0 with p+F>0 on the compact cells.
local averkov_negative_part_bound := proc(p, cells, x, {max_steps::posint := 200})
local F := 1, i;
    for i from 1 to max_steps while not univariate_positive_on_cells(p+F,cells,x,true) do
        F := 2*F;
    end do;
    return F;
end proc;

# Exponent n with T_n(1+2*eps/gamma)^2 >= K = (F+mu/2)/(2*eps*lambda).
# Floating point with a safety margin of one: the value only caps the
# search, and every accepted exponent is checked exactly.
local averkov_chebyshev_exponent := proc(gam, ep, lam, mu, F)
local K := (F+mu/2)/(2*ep*lam);
    if K<=1 then return 1; end if;
    return max(1, ceil(evalf[30](arccosh(sqrt(K))/arccosh(1+2*ep/gam)))+1);
end proc;

# cells intersected with {gens[j]+offsets[j]>=0} for every j.
local averkov_offset_cells := proc(cells0, gens, offsets, x)
local cells := cells0, j;
    for j from 1 to nops(gens) do
        if cells=[] then break; end if;
        cells := univariate_intersect_cells(cells,
                    univariate_select_cells(univariate_sign_chart(gens[j]+offsets[j],x),"nonnegative"));
    end do;
    return cells;
end proc;

# Smallest bound found with g<=2*gamma on cells: double, then bisect; when
# gamma=1 already bounds g, halve instead (capped: g<=0 on B bounds every gamma).
local averkov_gamma_bound := proc(g, cells, x, steps::nonnegint)
local c := 1, lo, mid, i;
    while not univariate_positive_on_cells(2*c-g,cells,x,false) do c := 2*c; end do;
    if c>1 then
        lo := c/2;
        for i from 1 to steps do
            mid := (lo+c)/2;
            if univariate_positive_on_cells(2*mid-g,cells,x,false) then c := mid; else lo := mid; end if;
        end do;
    else
        for i from 1 to 32 while univariate_positive_on_cells(c-g,cells,x,false) do c := c/2; end do;
    end if;
    return c;
end proc;

# Return SOS multipliers H such that f-dot_product(H,gens)>0 on {B_poly>=0}.
# Every bound and every accepted exponent is checked by exact sign decisions.
# Time is bounded by the caller's timelimit; max_N bounds the exponent search.
#
# Univariate refinement of Averkov's Lemma 7 with Chebyshev multipliers.
# Let T = {B>=0, f<=0} (finitely many cells). Only generators J that are
# negative somewhere on T can help, and each j in J gets its own constants:
#   phi_j(t) = lambda_j*t*T_{n_j}((t-gamma_j)/gamma_j)^2,  g_j <= 2*gamma_j on B,
#   lambda_j = mu/(4*|J|*gamma_j),  mu <= f on G = {B>=0, g_j>=-2*eps_j for j in J},
# where the eps_j are chosen so that G and T are disjoint. Then
#   * on G: f >= mu and sum phi_j(g_j) <= |J|*2*gamma_j*lambda_j = mu/2;
#   * off G, some g_j <= -2*eps_j and f-g >= -F-mu/2+2*eps_j*lambda_j*T_{n_j}(1+2*eps_j/gamma_j)^2,
# which is >= 0 for the exponent of averkov_chebyshev_exponent, of order
# sqrt(gamma_j/eps_j). The search runs over a level k: generator j uses
# max(1,ceil(k*n_j/K)), K = max n_j, so k=K reproduces the certified n_j
# and k is the largest exponent used (the quantity max_N caps).
local averkov_lemma_7 := proc(x, f, gens, B_poly, N_guess,
                              {max_N::nonnegint := 256, gamma_bisection_steps::nonnegint := 4})
local cells, bad_cells, J, j, k, i, gJ, mJ, eps0, eps, valid, lo, hi, mid;
local gam, mu, F, lam, n_cert, K_cert, exponents, mult, samples, values;
local N_top := 1, N_bot := 0, residual, sample_ok, exact_ok, next_N, H, e;
local N_cur, N_lo, N_hi;
$ifdef LOG_TIME
    INIT_START_LOG_TIME("averkov_lemma_7")
    try
$endif
    cells := univariate_select_cells(univariate_sign_chart(B_poly,x),"nonnegative");
    if univariate_positive_on_cells(f,cells,x) then
        return map(g_i -> 0, gens);
    end if;
    if max_N=0 then error "Averkov search exhausted max_averkov_N=%1",max_N; end if;

    # 1. The points to repair, and the generators that are negative there.
    bad_cells := univariate_intersect_cells(cells,
                    univariate_select_cells(univariate_sign_chart(f,x),"nonpositive"));
    J := select(j -> univariate_intersect_cells(bad_cells,
                        univariate_select_cells(univariate_sign_chart(gens[j],x),"negative"))<>[],
                [seq(i,i=1..nops(gens))]);
    if J=[] then error "f must be strictly positive on the nonnegative set of the generators"; end if;
    # Drop generators the others make redundant on T, highest degree first:
    # the certificate degree is max_j (2*n_j+1)*deg(g_j).
    for j in sort(J,(a,b) -> evalb(degree(gens[a],x)>degree(gens[b],x))) do
        if nops(J)>1 and averkov_offset_cells(bad_cells,map(i->gens[i],remove(`=`,J,j)),
                                              [0$(nops(J)-1)],x)=[] then
            J := remove(`=`,J,j);
        end if;
    end do;
    gJ := map(i->gens[i],J); mJ := nops(J);

    # 2. Separation margins: a common eps from findEps, then each eps_j
    # enlarged with the others fixed (validity is monotone in each eps_j).
    eps0 := findEps(x,gJ,f,B_poly);
    eps := [eps0$mJ];
    for k from 1 to mJ do
        valid := t -> evalb(averkov_offset_cells(bad_cells,gJ,
                              subsop(k=2*t,map(s->2*s,eps)),x)=[]);
        lo := eps[k]; hi := 2*lo;
        for i from 1 to 64 while valid(hi) do lo, hi := hi, 2*hi; end do;
        if valid(hi) then
            lo := hi;
        else
            for i from 1 to 3 do
                mid := (lo+hi)/2;
                if valid(mid) then lo := mid; else hi := mid; end if;
            end do;
        end if;
        eps := subsop(k=lo,eps);
    end do;

    # 3. Per-generator range bounds, positivity margin, scales and exponents.
    gam := map(g -> averkov_gamma_bound(g,cells,x,gamma_bisection_steps),gJ);
    mu := averkov_positive_lower_bound(f,averkov_offset_cells(cells,gJ,map(s->2*s,eps),x),x);
    F := averkov_negative_part_bound(f,cells,x);
    lam := [seq(mu/(4*mJ*gam[k]),k=1..mJ)];
    n_cert := [seq(averkov_chebyshev_exponent(gam[k],eps[k],lam[k],mu,F),k=1..mJ)];
    K_cert := max(op(n_cert));
    DEBUG(__FILE__, __LINE__, lprint(">> Univariate Averkov J, gamma, eps, mu, F, n",J,gam,eps,mu,F,n_cert));

    exponents := level -> [seq(max(1,ceil(level*n_cert[k]/K_cert)),k=1..mJ)];
    # Multiplier of generator k at exponent n, at a value or polynomial w of g_k.
    mult := (k,n,w) -> lam[k]*averkov_chebyshev(n,(w-gam[k])/gam[k])^2;
    residual := proc(level)
    local ex := exponents(level), k2;
        return expand(f-add(mult(k2,ex[k2],gJ[k2])*gJ[k2],k2=1..mJ));
    end proc;
    # A nonpositive value at a point of {B>=0} rejects a level without
    # expanding the residual. Passing samples is necessary, not sufficient.
    samples := averkov_sample_points(cells,B_poly,x);
    values := map(s -> [subs(x=s,f), map(g -> subs(x=s,g),gJ)], samples);
    sample_ok := proc(level)
    local ex := exponents(level), v, k2;
        for v in values do
            if v[1]-add(mult(k2,ex[k2],v[2][k2])*v[2][k2],k2=1..mJ) <= 0 then return false; end if;
        end do;
        return true;
    end proc;
    next_N := proc(n)
    local m;
        if n>=max_N then
            error "Averkov search exhausted max_averkov_N=%1",max_N;
        end if;
        m := min(2*n,max_N);
        # Never jump past the certified level; beyond it (only if rounding
        # made the bound too small) plain doubling resumes.
        if n<K_cert then m := min(m,K_cert); end if;
        # Retain the configured guess as an extra probe within the cap.
        if n<N_guess and N_guess<m then m := N_guess; end if;
        return m;
    end proc;
    exact_ok := n -> univariate_positive_on_cells(residual(n),cells,x);
    # Keep N_bot rejected and N_top sample-admissible; resume doubling when
    # N_top fails the exact check. No monotonicity assumption is needed for
    # correctness, and the returned level need not be the smallest one.
    do
        while not sample_ok(N_top) do
            N_bot := N_top;
            N_top := next_N(N_top);
        end do;
$ifdef ENABLE_BINARY_SEARCH_AVKL
        # Bisect (N_bot,N_top] by samples alone. When samples are informative,
        # one exact check accepts the smallest sample-admissible level.
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
            # rejected N_lo and the verified N_top with exact checks.
            while N_top-N_lo>1 do
                N_cur := iquo(N_top+N_lo,2);
                if sample_ok(N_cur) and exact_ok(N_cur) then N_top := N_cur; else N_lo := N_cur; end if;
            end do;
$endif
            break;
        end if;
        N_bot := N_top;
        N_top := next_N(N_top);
    end do;
    e := exponents(N_top);
    DEBUG(__FILE__, __LINE__, lprint(">> Univariate Averkov level, exponents",N_top,e));
    H := map(g_i -> 0, gens);
    for k from 1 to mJ do
        H := subsop(J[k]=mult(k,e[k],gJ[k]),H);
    end do;
    return H;
$ifdef LOG_TIME
    finally
        END_LOG_TIME("averkov_lemma_7");
    end try;
$endif
end proc;
