# Generates tests/benchmarks/lbp_stress.mpl (run from the repository root).
#
# Each case [f, basis, arch] satisfies, checked exactly:
# - basis has at most 4 generators of degree <= 3 and contains the
#   Archimedean generator arch = R^2 - x^2, so the quadratic module is
#   Archimedean and every f > 0 on S = {basis >= 0} is a member (Putinar);
# - S has nonempty interior;
# - deg f in 3..7 exceeds deg(arch) = 2 and f > 0 on S, with f(x) >= delta
#   up to a rational rounding of the minimum, delta in {1, 1/10, 1/100};
# - f < 0 somewhere on R, so spCertificates cannot return early.
interface(echo=0, warnlevel=0):
with(RootFinding, Isolate):
randomize(20261113):
NUM_CASES := 30:
rand_point := proc(R) rand(-2*R+1 .. 2*R-1)()/2 end proc:

random_generator := proc(R)
local a, b, c, kind := rand(1..4)();
    a, b, c := op(sort([rand_point(R), rand_point(R), rand_point(R)]));
    if kind = 1 then
        return `if`(rand(0..1)() = 0, x - a, a - x);
    elif kind = 2 then
        return (x - a)*(x - c);          # removes (a, c)
    elif kind = 3 then
        return -(x - a)*(x - c);         # keeps [a, c]
    else
        return `if`(rand(0..1)() = 0, 1, -1)*(x - a)*(x - b)*(x - c);
    end if;
end proc:

min_over_S := proc(p, basis)
local candidates, inS;
    candidates := map(rhs, Isolate(diff(p, x), x, digits=30));
    candidates := [op(candidates), op(map(q -> op(map(rhs, Isolate(q, x, digits=30))), basis))];
    inS := select(t -> andmap(q -> evalf(subs(x=t, q)) >= -1e-20, basis), candidates);
    return min(op(map(t -> evalf(subs(x=t, p)), inS)));
end proc:

cases := []:
k := 0:
while nops(cases) < NUM_CASES do
    R := rand(2..5)();
    arch := R^2 - x^2;
    basis := [arch, seq(random_generator(R), i = 1 .. rand(1..3)())];
    basis := combinat:-randperm(basis);
    S := SolveTools:-SemiAlgebraic(map(q -> q > 0, basis), [x]);
    if S = [] then next; end if;
    d := rand(3..7)();
    p := randpoly(x, degree = d, coeffs = rand(-5..5));
    if degree(p, x) <> d then next; end if;
    delta := [1, 1/10, 1/100][irem(nops(cases), 3) + 1];
    f := expand(p - floor(1000*min_over_S(p, basis))/1000 + delta);
    if SolveTools:-SemiAlgebraic([f <= 0, op(map(q -> q >= 0, basis))], [x]) <> [] or
       SolveTools:-SemiAlgebraic([f < 0], [x]) = [] then
        next;
    end if;
    cases := [op(cases), [f, basis, arch]];
end do:
external_benchmarks := cases:
save external_benchmarks, "tests/benchmarks/lbp_stress.mpl":
printf("Saved %d cases; degrees of f: %a; basis sizes: %a\n", nops(cases),
       map(c -> degree(c[1], x), cases), map(c -> nops(c[2]), cases)):
quit:
