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
basis using the SOS representation of `B`. This replaces the calls to
`lower_bound_poly` and `averkov_extended_lemma` in `spCertificates`.
The return format remains `[sigma0, sigma1, ..., sigman]`, satisfying
`f = sigma0 + add(sigmai*basis[i], i=1..nops(basis))`.

For a single compact generator, the monogenic method is also available directly:

    monogenicCertificates(10^9*(1-x^2)+1, 1-x^2, x);
    # [1, 1000000000]

Both `spCertificates` and `monogenicCertificates` accept the keyword limits
`max_k=100` and `max_bisection=200`. Exhausted limits raise an error;
they do not indicate that no certificate exists. Inputs to the monogenic
method must be univariate polynomials with rational coefficients, with
`g` having even degree and a negative leading coefficient, and `f` strictly
positive on `{g >= 0}`. Empty and singleton feasible sets are supported.

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
