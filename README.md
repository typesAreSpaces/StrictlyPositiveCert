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

After constructing a compact bound polynomial `B` and applying
`averkov_lemma_7`, `spCertificates` now writes the remaining polynomial as
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
budget. These limits are cooperative Maple time limits. Explicit SOS square
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
