# Exact epsilon separation and Averkov regressions; read source, not an archive.
read "src/StrictlyPositiveCert.mpl":
with(SolveTools, SemiAlgebraic):
kernelopts(opaquemodules=false):

check_epsilon := proc()
local x, c := 73205080757/10^11, item, f, gens, eps, H, h, i, caught, original_monogenic, original_positive,
      original_samples, B, exact_checks, N;
    for item in [
        [x-3/4,[x-1,2*x-1,3*x-1,4-x^2]],
        [-1,[x-c,2-x^2-x-c]],
        [x+1/10^10,[1-x^2,x]],
        [x+1,[-x^2]],
        [1,[1-x^2]],
        [x^2-1/2,[-(x^2-1)*(x^2-4)]]
    ] do
        f,gens := op(item);
        eps := StrictlyPositiveCert:-findEps(x,gens,f);
        if not type(eps,rational) or eps<=0 or
        SemiAlgebraic([f<=0,op(map(g->g+2*eps>=0,gens))],[x])<>[] then
            error "unsafe epsilon for %1: %2",item,eps;
        end if;
    end do;
    caught := false;
    try StrictlyPositiveCert:-findEps(x,[1-x^2],x);
    catch "f must be strictly positive": caught := true;
    end try;
    if not caught then error "epsilon search accepted invalid positivity input"; end if;

    # Exact positivity on B must not use rounded algebraic endpoints.
    # Small positive margins must not enter floating logarithms or solve.
    for item in [
        [14142135623/10^10-x,[2-x^2,1-x]],
        [x+1/10^10,[1-x^2,x]],
        [x+1/2,[1-x^2,x]],
        [-1,[-x^2,x^2-1]]
    ] do
        f,gens := op(item);
        H := timelimit(20,StrictlyPositiveCert:-spCertificates(f,gens,x));
        if nops(H)<>nops(gens)+1 or
        expand(f-StrictlyPositiveCert:-dot_product([1,op(gens)],H))<>0 then
            error "certificate identity failed: %1",item;
        end if;
        for h in H do
            if not type(h,polynom(rational,x)) or SemiAlgebraic([h<0],[x])<>[] then
                error "invalid SOS multiplier: %1",item;
            end if;
        end do;
    end do;
    # The gamma bound is checked on B, rather than globally: high even-degree
    # positive generators can be unbounded on R while bounded on the compact set.
    H := StrictlyPositiveCert:-averkov_lemma_7(x,x+1/2,[1-x^2,x,x^4],1-x^2,70);
    if SemiAlgebraic([1-x^2>=0,
        x+1/2-StrictlyPositiveCert:-dot_product([1-x^2,x,x^4],H)<=0],[x])<>[] then
        error "Averkov residual is not strictly positive on B";
    end if;
    # eps* = 1/4 is itself invalid: halving reaches 1/8, then three
    # bisection steps keep the largest valid midpoint below eps*.
    eps := StrictlyPositiveCert:-findEps(x,[1-x^2,x],x+1/2);
    if eps <> 15/64 then error "bisection did not refine eps: %1",eps; end if;
    # eps* = 101/2: doubling stops at 64, bisection keeps 48.
    eps := StrictlyPositiveCert:-findEps(x,[-x^2-100],x+1);
    if eps <> 48 then error "doubling did not grow eps: %1",eps; end if;
    if StrictlyPositiveCert:-findEps(x,[1-x^2],1) <> 1 then
        error "empty separation set did not return immediately";
    end if;

    # Left-benchmark instance: a near-unit eps made the exponent search time out.
    f := x+241/50;
    gens := [(x-112/25)*(x-294/25)*(x-343/25), (-x+112/25)*(x-294/25)*(x-343/25), 154449/625-x^2];
    H := timelimit(60,StrictlyPositiveCert:-spCertificates(f,gens,x));
    if expand(f-StrictlyPositiveCert:-dot_product([1,op(gens)],H))<>0 or
    ormap(h -> SemiAlgebraic([h<0],[x])<>[], H) then
        error "left benchmark certificate is invalid";
    end if;

    # A residual already nonnegative on R must not reach the monogenic stage.
    original_monogenic := eval(StrictlyPositiveCert:-monogenicCertificates);
    unprotect(StrictlyPositiveCert:-monogenicCertificates);
    try
        StrictlyPositiveCert:-monogenicCertificates := proc()
            error "monogenic stage reached";
        end proc;
        f := (x-53/50)*(-x+53/50)+1;
        gens := [x-53/50, -x+53/50, 23409/2500-x^2];
        H := StrictlyPositiveCert:-spCertificates(f,gens,x);
    finally
        StrictlyPositiveCert:-monogenicCertificates := eval(original_monogenic);
        protect(StrictlyPositiveCert:-monogenicCertificates);
    end try;
    if expand(f-StrictlyPositiveCert:-dot_product([1,op(gens)],H))<>0 or
    SemiAlgebraic([H[1]<0],[x])<>[] then
        error "nonnegative residual shortcut produced an invalid certificate";
    end if;
    # Uninformative samples must not reduce the exponent search to one
    # exact check per exponent: the count must stay below the exponent found.
    f := x+241/50;
    gens := [(x-112/25)*(x-294/25)*(x-343/25), (-x+112/25)*(x-294/25)*(x-343/25), 154449/625-x^2];
    B := StrictlyPositiveCert:-bound_poly(gens,x)[1];
    original_positive := eval(StrictlyPositiveCert:-univariate_positive_on_cells);
    original_samples := eval(StrictlyPositiveCert:-averkov_sample_points);
    unprotect(StrictlyPositiveCert:-univariate_positive_on_cells);
    unprotect(StrictlyPositiveCert:-averkov_sample_points);
    exact_checks := 0;
    try
        StrictlyPositiveCert:-averkov_sample_points := proc() [] end proc;
        # Residual checks are the three-argument calls on nonlinear polynomials.
        StrictlyPositiveCert:-univariate_positive_on_cells := proc(p, cells, y)
            if nargs=3 and degree(expand(p),y)>1 then exact_checks := exact_checks+1; end if;
            return original_positive(args);
        end proc;
        H := StrictlyPositiveCert:-averkov_lemma_7(x,f,gens,B,70);
    finally
        StrictlyPositiveCert:-univariate_positive_on_cells := eval(original_positive);
        StrictlyPositiveCert:-averkov_sample_points := eval(original_samples);
        protect(StrictlyPositiveCert:-univariate_positive_on_cells);
        protect(StrictlyPositiveCert:-averkov_sample_points);
    end try;
    N := degree(H[3],x)/(2*degree(gens[3],x));
    if SolveTools:-SemiAlgebraic([B>=0,f-StrictlyPositiveCert:-dot_product(H,gens)<=0],[x])<>[] then
        error "Averkov residual without samples is not strictly positive on B";
    end if;
    if exact_checks>=N then
        error "exponent search used %1 exact checks for exponent %2",exact_checks,N;
    end if;
    printf("PASS: exact dyadic epsilon, independent feasibility oracle, invalid inputs, and Averkov precision regressions\n");
end proc:
check_epsilon():
