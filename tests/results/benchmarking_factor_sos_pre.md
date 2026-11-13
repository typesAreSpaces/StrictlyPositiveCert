# FACTOR_SOS_PRE comparison on tests/benchmarking.mpl

The benchmark file and its helpers were run unchanged. The only active suite is
`strictly_positive_polynomials_stengle_benchmark`: f=1-x^2+d, basis=[(1-x^2)^k],
d in [10, 5, 1, 1/2, 1/3], k in [1, 3]. Its existing 300-second limits remain.
Both configurations use the current Makefile settings: Averkov binary search and
debugging enabled, EPS_FACTOR=17/20, N_GUESS_AVKL=70, N_GUESS_LS=600.
WEIFENG_OPTIMIZATION and other optional flags remain disabled. The helper calls
the package method “Weifeng” regardless of this optimization flag.

Maple 2022 loaded package source through `tests/run_factor_sos_benchmarking.mpl`
before reading the actual benchmark file from the tests directory. No installed
.mla, benchmark case, or helper was changed. Debug output was redirected in both
configurations. Seven fresh-process runs per configuration were executed
sequentially, alternating order after the initial off/on pair. Tables use median
CPU seconds reported by the original helpers (not wall time); loading and output
after each timed construction are excluded. This differs from the previous
small-test harness, which measured warmed wall times. Degrees are the maximum
coefficient degrees reported by the existing `quickdegree` helper, not degrees
of individual square summands. Neither original helper independently validates
the certificate identity or coefficient nonnegativity during this benchmark.

## StrictlyPositiveCert results

Negative changes mean faster with factoring enabled.

| d | k | Off (s) | On (s) | Difference (s) | Change | Degree off → on | CPU range off / on (s) |
|---|---:|---:|---:|---:|---:|---:|---|
| 10 | 1 | 0.125 | 0.127 | +0.002 | +1.6% | 2 → 2 | 0.092–0.167 / 0.114–0.158 |
| 10 | 3 | 0.037 | 0.043 | +0.006 | +16.2% | 6 → 6 | 0.033–0.039 / 0.029–0.054 |
| 5 | 1 | 0.028 | 0.030 | +0.002 | +7.1% | 2 → 2 | 0.027–0.038 / 0.026–0.036 |
| 5 | 3 | 0.033 | 0.038 | +0.005 | +15.2% | 6 → 6 | 0.026–0.042 / 0.029–0.042 |
| 1 | 1 | 0.034 | 0.029 | -0.005 | -14.7% | 2 → 2 | 0.027–0.042 / 0.024–0.042 |
| 1 | 3 | 0.036 | 0.040 | +0.004 | +11.1% | 6 → 6 | 0.032–0.043 / 0.028–0.049 |
| 1/2 | 1 | 0.031 | 0.035 | +0.004 | +12.9% | 2 → 2 | 0.026–0.033 / 0.018–0.039 |
| 1/2 | 3 | 0.031 | 0.053 | +0.022 | +71.0% | 6 → 6 | 0.029–0.036 / 0.042–0.084 |
| 1/3 | 1 | 0.041 | 0.028 | -0.013 | -31.7% | 2 → 2 | 0.033–0.047 / 0.020–0.034 |
| 1/3 | 3 | 0.371 | 0.403 | +0.032 | +8.6% | 1386 → 1386 | 0.334–0.405 / 0.361–0.417 |

Sum of per-case medians: **0.767 → 0.826 s**, **+7.7%**. This is not total process wall time.

## RealCertify results

RealCertify is a separate package and does not use FACTOR_SOS_PRE. Its timings
provide context, not evidence of an effect of this flag. The benchmark prints
“Succeeds RealCertify” even when the helper returns no certificate; the status
below is taken from the helper's actual output.

| d | k | Off median (s) | On median (s) | Degree off → on | Status off / on |
|---|---:|---:|---:|---:|---|
| 10 | 1 | 0.073 | 0.070 | 4 → 4 | success / success |
| 10 | 3 | 0.009 | 0.009 | 6 → 6 | success / success |
| 5 | 1 | 0.009 | 0.010 | 4 → 4 | success / success |
| 5 | 3 | 0.010 | 0.011 | 6 → 6 | success / success |
| 1 | 1 | 0.008 | 0.009 | 4 → 4 | success / success |
| 1 | 3 | 0.010 | 0.010 | 6 → 6 | success / success |
| 1/2 | 1 | 0.010 | 0.010 | 4 → 4 | success / success |
| 1/2 | 3 | 0.011 | 0.010 | 6 → 6 | success / success |
| 1/3 | 1 | 0.049 | 0.009 | 4 → 4 | success / success |
| 1/3 | 3 | — | — | — → — | no_certificate / no_certificate |

Small timing differences are noisy; min/max ranges are included above. Raw
per-run results are in `benchmarking_factor_sos_pre_runs.json`; logs remain in
`/tmp/factor_sos_benchmarking_{off,on}*.log`.

## Reproduce

Run from the repository root, adding only `-DFACTOR_SOS_PRE` for the enabled run:

```sh
maple -q -DENABLE_BINARY_SEARCH_AVKL -DENABLE_DEBUGGING -DEPS_FACTOR=17/20 -DN_GUESS_AVKL=70 -DN_GUESS_LS=600 tests/run_factor_sos_benchmarking.mpl
maple -q -DENABLE_BINARY_SEARCH_AVKL -DENABLE_DEBUGGING -DEPS_FACTOR=17/20 -DN_GUESS_AVKL=70 -DN_GUESS_LS=600 -DFACTOR_SOS_PRE tests/run_factor_sos_benchmarking.mpl
```
