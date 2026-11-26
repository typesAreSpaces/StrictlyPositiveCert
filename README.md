# StrictlyPositiveCert

The package `StrictlyPositiveCert` computes certificates of strictly positive
univariate polynomials over basic semialgebraic sets with an Archimedean
quadratic module. Its construction is based on `Algorithm 4` in (Shang et
al. 2025), with an exact monogenic certificate search as the final stage.

## Installation

Use the `make` command at the root of the project to produce the maple
library file `StrictlyPositiveCert.mla`. Then, load the
`StrictlyPositiveCert` package by
adding the following two lines at the beginning of your maple file:

    libname := libname, "<PATH-TO-STRICTLYPOSITIVECERT>/StrictlyPositiveCert.mla":
    with(StrictlyPositiveCert):

When built with `-DENABLE_DEBUGGING` (for example,
`make DEBUG_MODE=-DENABLE_DEBUGGING`), the library prints its source branch
when loaded with `with(StrictlyPositiveCert)`. The branch name and debug
setting are saved in the `.mla`; rebuild the library to change them.

## Usage and examples

The command `spCertificate` has the following inputs:

- `f`: Denotes the input polynomial
- `basis`: This is a list of polynomials (generators) in the free
  variable `x`
- `x`: Denotes the free variable of the ring of univariate polynomial

The directory `tests` contains benchmarks and examples reported in
(Shang et al. 2025).

For a simple example, we can execute

    libname := libname, "../StrictlyPositiveCert.mla":
    with(StrictlyPositiveCert);

    spCertificates(x+1, [-x^2], x);

The output is


    > libname := libname, "../StrictlyPositiveCert.mla":
    > with(StrictlyPositiveCert);
                                                            [bound_info, dot_product, findEps, spCertificates]


    > spCertificates(x+1, [-x^2], x);
                                                                                         2
                                                                          [x + 1 + 9/10 x , 9/10]

We can verify $$x+1 = (x+1+\frac{9}{10}x^2 + \frac{9}{10} \dot (-x^2))$$ and
$$x+1+\frac{9}{10}x^2 = (\frac{1}{2}x+1)^2+\frac{13}{20}x^2$$ is a sums of squares.

## References



Shang, Weifeng, Chenqi Mou, Jose Abel Castellanos Joo, and Deepak Kapur. 2025. *Computing Certificates of Strictly Positive Polynomials in
Archimedean Quadratic Modules*. <https://arxiv.org/abs/2503.11119>.



## Monogenic certificate stage

After constructing a compact bound polynomial `B`, `spCertificates` first tries
exact univariate reductions and then uses a bounded `averkov_lemma_7` fallback.
It writes the remaining polynomial as
`f2 = sigma0 + sigma1*B` using the exact search ported from
`monogenicCert.wl`. It searches over

    sigma1 = lambda*((M-B)/M)^(2*k)

with positive rational `lambda`. Real-root counting and exact semialgebraic
feasibility checks verify that `sigma0 = f2-sigma1*B` is strictly positive
on the real line. The multiplier is then distributed over the original
basis using the SOS representation of `B`.
The return format remains `[sigma0, sigma1, ..., sigman]`, satisfying
`f = sigma0 + add(sigmai*basis[i], i=1..nops(basis))`.

For a single compact generator, the monogenic method is also available directly:

    monogenicCertificates(10^9*(1-x^2)+1, 1-x^2, x);
    # Returns [sigma0, sigma1]; the first verified search returns its certificate.

Both `spCertificates` and `monogenicCertificates` accept `max_k=100` and
`max_bisection=200`. The search advances every active multiplier search once per round. Exponent
`k` is admitted by round `k^2` (earlier if no searches remain active), limiting
high-degree work while lower-degree searches are still progressing. Each exponent retains its own
multiplier and brackets; it does not have to exhaust its budget before the
next exponent is tried. `max_bisection` retains the initial trial plus that
many updates for each exponent. The first verified certificate is returned;
interleaving and timeouts can change its degree compared with a depth-first
search. The output is not guaranteed to have the smallest degree.

A trial starts with `query_time_limit=1/10` second. On a timeout, the same trial
is deferred to the next round with twice its previous allowance, up to
`max_query_retries=4` retries. A timeout never establishes infeasibility or
moves the multiplier bracket. `max_search_time=300` bounds elapsed search
time, including normalization; each trial is also capped by the remaining
budget. For `spCertificates`, one deadline is shared across preprocessing,
univariate reduction, Averkov, and the final monogenic search; the latter gets
only the remaining allowance. These limits are cooperative Maple time limits. Explicit SOS square
decomposition has its own iteration budget and is outside the search deadline.
Exhausted budgets raise an error; they do not indicate that no certificate
exists. Inputs to the monogenic method must be univariate polynomials with
rational coefficients, with `g` having even degree and a negative leading
coefficient, and `f` strictly positive on `{g >= 0}`. Empty and singleton
feasible sets are supported.

The default output contains polynomial multipliers, which are SOS because
strictly positive univariate real polynomials are SOS. To obtain explicit
weighted rational squares, use:

    monogenicCertificates(x+3/2, 1-x^2, x, output_sos=true,
                          max_sos_iterations=66);

In this mode each multiplier is represented by a list of `[c, q]` pairs
meaning `add(c*q^2)`. The optional square decomposition uses Schweighofer's
quadratic reduction, as in the supplied Wolfram implementation. Its separate
`max_sos_iterations` limit defaults to 64 and applies to each reduction step.
Increasing it may be necessary for large coefficients. The default polynomial
output does not require this decomposition and is not subject to that limit.

Run the source regressions with `make test` (Maple and Python 3 required).
The tests check polynomial identities, nonnegativity, explicit square
representations, search limits, and the public multi-generator integration.

## Exact univariate sign decisions

The certificate routines use exact real-root isolation and Sturm counts in place
of general `SemiAlgebraic` calls. The monogenic search constructs the sign cells
of its fixed generator once and reuses them for normalization and lambda trials.
Its dyadic bounds retain the same equality conventions and minimum cap as before.

Strict positivity rejects every included zero. Nonnegativity instead uses the
odd-multiplicity part from a square-free decomposition, allowing even-multiplicity
zeros. Sign cells retain exact algebraic endpoints and distinguish open intervals,
closed intervals, isolated points, and unbounded components. Averkov and epsilon
computations obtain their interval sets from these same sign constructions.

The regression suite compares set boundaries, positivity decisions, lambda
classification, and dyadic bounds against `SemiAlgebraic` as an independent oracle.

## Exact Averkov parameter search

The Averkov stage checks strict positivity on its compact bound using exact
sign cells and Sturm counts. It chooses a positive dyadic `gamma` such that
all generators are at most `2*gamma` on that set. Its epsilon search looks for
a large valid `eps`, meaning no point satisfies `f <= 0` and every
`g_i + 2*eps >= 0`; the original feasible set is checked first to reject invalid
positivity inputs. Validity is monotone, so `findEps` doubles from 1 while valid
(at most 64 doublings) or halves until valid, then runs three bisection steps
between the last valid and first invalid values and returns the last valid one.
It returns 1 immediately when there is nothing to separate. A larger `eps`
makes the Averkov exponent converge much faster.

The prefactor is fixed at one. Exponents increase, up to `max_averkov_N=256`.
Candidates are first rejected cheaply by evaluating the residual on a rational
grid of points in `{B >= 0}`; doubling and the optional binary search use only
this test. The final exponent is then verified strictly positive on the bound
with exact sign decisions, and doubling resumes if that check fails. No
monotonicity is assumed.
No floating-point minimization, maximization, logarithms, or symbolic alpha
optimization is used in this stage. Certificate degrees and coefficient sizes
may differ, and the exact search does not promise a minimum-degree certificate.
`N_GUESS_AVKL` remains an optional extra probe
when exponent doubling crosses the configured guess.

## Symbolic univariate reductions

Before the Averkov fallback, `spCertificates` tries polynomial division by each
nonconstant generator. A quotient nonnegative on the whole real line is an SOS
multiplier, including quotients such as `x^2` with zeros. The remainder need only
be strictly positive on `{B >= 0}`. This is an intermediate reduction, not a
change to the exported `quickCertificates` contract. When the residual is already
globally nonnegative, the complete certificate is returned immediately; otherwise
it proceeds to the monogenic stage.

The next heuristic searches for a positive rational `lambda` in

    f - lambda*P > 0 on {B >= 0},  P = add(weights[i]*gens[i]).

Weights are fixed squares chosen from `1`, `x^2`, `(x-1)^2`, and `(x+1)^2`,
with joint directions and one-generator directions tried in that order. A point
with `P >= 0` and `f <= 0`, including a zero of `P`, rejects the direction
immediately. Exact sign decisions on
`{P > 0}` and `{P < 0}` classify a trial as too large, too small, infeasible, or
verified. Brackets and all scales are rational; no numerical optimization or SDP
solver is used. Coefficient-size normalization initializes the scale so tiny
constant directions can be useful without a long doubling search.

New `spCertificates` options are:

- `use_univariate_heuristics=true`: enable division and scalar reductions.
- `max_heuristic_trials=64`: total scalar trial budget across directions; zero
  skips scalar searches while retaining division.
- `max_scalar_updates=32`: at most this many bracket updates plus the initial
  trial for each direction, subject to the total trial budget.
- `max_averkov_N=256`: cap the fallback exponent before polynomial construction.
  Zero allows reductions and the zero-multiplier shortcut, but disables positive
  Averkov exponents.

Heuristics may miss a certificate and then fall back to Averkov. Exhausting the
shared deadline or exponent cap raises a stage-specific error, not a claim of
infeasibility. `max_k`, `max_bisection`, `query_time_limit`, and `max_query_retries`
continue to control only the monogenic stage. Constant filtering preserves the
original generator indices without restarting the deadline. Inputs must be
rational univariate polynomials and generators must be supplied as a list.

For example:

    spCertificates(x^3+1/10, [1-x^2,x], x, max_averkov_N=0);
    # [1/10, 0, x^2]
    spCertificates(-1, [-x^2,x^2-1/10^6], x, max_averkov_N=0);
    # [1, 2000000, 2000000]

These limits remain cooperative Maple limits; a native operation can overrun its
allowance. The fixed square directions and the exponent cap do not guarantee
minimum-degree output. Standalone `monogenicCertificates` behavior is unchanged.

## Univariate Chebyshev construction in `averkov_lemma_7` (experiment_univheuristics)

`averkov_lemma_7` uses a univariate refinement of Averkov's Lemma 7 with
Chebyshev multipliers `lambda_j*T_{n_j}((g_j-gamma_j)/gamma_j)^2`:

- Only generators that are negative somewhere on `{B>=0, f<=0}` are used.
  Generators made redundant there by the others are dropped, highest degree
  first, and get multiplier 0.
- Each remaining generator has its own range bound `gamma_j`, separation
  margin `eps_j`, scale `lambda_j = mu/(4*|J|*gamma_j)` and certified exponent
  `n_j ~ sqrt(gamma_j/eps_j)`. The margins start from `findEps` and are
  enlarged one at a time.
- The search runs over a level `k`. Generator `j` uses `max(1, ceil(k*n_j/K))`
  with `K = max n_j`, and `max_averkov_N` caps `k`, which is the largest
  exponent used.

All constants come from exact univariate sign decisions, and every accepted
level is verified exactly.
