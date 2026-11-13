# Factoring and stage division: four-configuration benchmark

Ran the unchanged active suite in `tests/benchmarking.mpl`: ten Stengle cases
f=1-x^2+d, basis=[(1-x^2)^k], d in [10, 5, 1, 1/2, 1/3], k in [1, 3].
The file's existing 300-second limits and RealCertify calls remain unchanged.
Only FACTOR_SOS_PRE and ENABLE_STAGE_DIVISION vary. All other compiler flags
and values come from the current Makefile's dry-run command:

```text
maple  -DENABLE_BINARY_SEARCH_AVKL  -DENABLE_DEBUGGING  -DEPS_FACTOR=17/20 -DN_GUESS_AVKL=70 -DN_GUESS_LS=600 ./scripts/installation_script.mpl
```

Maple 2022; seven fresh processes per combination, run sequentially with
rotating configuration order. Source is loaded through the existing wrapper,
then the actual benchmark file is read from the tests directory. Installed
.mla files and benchmark source/helpers are untouched. Source hashes were
checked throughout; they and per-run results are saved in the JSON artifact.
Debug output remains enabled and is redirected to logs for every run.

Times are median CPU seconds reported by the existing checkWeifeng helper,
excluding source loading and subsequent certificate printing. Degrees are
the maximum certificate coefficient degrees reported by the existing
quickdegree helper, not degrees of individual square summands. No extra
certificate verification flags were enabled. The original benchmark does
not independently verify identities/nonnegativity. Small timing differences
are noisy; per-case ranges appear below.

Columns identify the two flags as F=FACTOR_SOS_PRE, D=ENABLE_STAGE_DIVISION.

## StrictlyPositiveCert median CPU times (seconds)

| d | k | F0 D0 | F1 D0 | F0 D1 | F1 D1 |
|---|---:|---:|---:|---:|---:|
| 10 | 1 | 0.126 | 0.123 | 0.115 | 0.111 |
| 10 | 3 | 0.031 | 0.037 | 0.054 | 0.056 |
| 5 | 1 | 0.030 | 0.031 | 0.031 | 0.028 |
| 5 | 3 | 0.034 | 0.036 | 0.035 | 0.043 |
| 1 | 1 | 0.025 | 0.036 | 0.029 | 0.032 |
| 1 | 3 | 0.035 | 0.040 | 0.045 | 0.046 |
| 1/2 | 1 | 0.029 | 0.030 | 0.028 | 0.031 |
| 1/2 | 3 | 0.037 | 0.078 | 0.086 | 0.093 |
| 1/3 | 1 | 0.036 | 0.031 | 0.028 | 0.029 |
| 1/3 | 3 | 0.360 | 0.374 | 0.379 | 0.388 |

## Certificate coefficient degrees

| d | k | F0 D0 | F1 D0 | F0 D1 | F1 D1 |
|---|---:|---:|---:|---:|---:|
| 10 | 1 | 2 | 2 | 0 | 0 |
| 10 | 3 | 6 | 6 | 6 | 6 |
| 5 | 1 | 2 | 2 | 0 | 0 |
| 5 | 3 | 6 | 6 | 6 | 6 |
| 1 | 1 | 2 | 2 | 0 | 0 |
| 1 | 3 | 6 | 6 | 6 | 6 |
| 1/2 | 1 | 2 | 2 | 0 | 0 |
| 1/2 | 3 | 6 | 6 | 6 | 6 |
| 1/3 | 1 | 2 | 2 | 0 | 0 |
| 1/3 | 3 | 1386 | 1386 | 1386 | 1386 |

## Aggregate

| Configuration | Sum of per-case medians (s) | Change vs F0 D0 |
|---|---:|---:|
| F0 D0 | 0.743 | +0.0% |
| F1 D0 | 0.816 | +9.8% |
| F0 D1 | 0.830 | +11.7% |
| F1 D1 | 0.857 | +15.3% |

These aggregates are sums of per-case medians, not total process wall times.

Package failures: none; all 280 package calls completed successfully.

## Package CPU ranges (minimum–maximum seconds)

| d | k | F0 D0 | F1 D0 | F0 D1 | F1 D1 |
|---|---:|---:|---:|---:|---:|
| 10 | 1 | 0.114–0.151 | 0.118–0.145 | 0.101–0.132 | 0.082–0.135 |
| 10 | 3 | 0.025–0.037 | 0.021–0.042 | 0.044–0.067 | 0.047–0.065 |
| 5 | 1 | 0.028–0.038 | 0.028–0.041 | 0.020–0.035 | 0.026–0.032 |
| 5 | 3 | 0.029–0.040 | 0.019–0.039 | 0.021–0.045 | 0.031–0.063 |
| 1 | 1 | 0.024–0.038 | 0.017–0.041 | 0.016–0.040 | 0.021–0.036 |
| 1 | 3 | 0.029–0.042 | 0.033–0.049 | 0.019–0.055 | 0.044–0.053 |
| 1/2 | 1 | 0.020–0.037 | 0.028–0.034 | 0.023–0.037 | 0.026–0.036 |
| 1/2 | 3 | 0.026–0.046 | 0.058–0.088 | 0.084–0.095 | 0.083–0.102 |
| 1/3 | 1 | 0.033–0.045 | 0.024–0.040 | 0.025–0.036 | 0.014–0.034 |
| 1/3 | 3 | 0.314–0.393 | 0.302–0.398 | 0.331–0.392 | 0.339–0.417 |

## RealCertify

This separate package does not use either flag. Each cell below gives median
CPU seconds and certificate degree for its successful calls, or failure status.
The script's unconditional “Succeeds RealCertify” messages are not treated as
proof of success; the parser checks for the helper's no-certificate output.

| d | k | F0 D0 | F1 D0 | F0 D1 | F1 D1 |
|---|---:|---|---|---|---|
| 10 | 1 | 0.070 s, degree 4 (success) | 0.069 s, degree 4 (success) | 0.067 s, degree 4 (success) | 0.071 s, degree 4 (success) |
| 10 | 3 | 0.009 s, degree 6 (success) | 0.009 s, degree 6 (success) | 0.009 s, degree 6 (success) | 0.010 s, degree 6 (success) |
| 5 | 1 | 0.009 s, degree 4 (success) | 0.009 s, degree 4 (success) | 0.010 s, degree 4 (success) | 0.008 s, degree 4 (success) |
| 5 | 3 | 0.010 s, degree 6 (success) | 0.010 s, degree 6 (success) | 0.010 s, degree 6 (success) | 0.009 s, degree 6 (success) |
| 1 | 1 | 0.010 s, degree 4 (success) | 0.009 s, degree 4 (success) | 0.010 s, degree 4 (success) | 0.008 s, degree 4 (success) |
| 1 | 3 | 0.010 s, degree 6 (success) | 0.009 s, degree 6 (success) | 0.010 s, degree 6 (success) | 0.009 s, degree 6 (success) |
| 1/2 | 1 | 0.010 s, degree 4 (success) | 0.011 s, degree 4 (success) | 0.011 s, degree 4 (success) | 0.009 s, degree 4 (success) |
| 1/2 | 3 | 0.012 s, degree 6 (success) | 0.010 s, degree 6 (success) | 0.010 s, degree 6 (success) | 0.008 s, degree 6 (success) |
| 1/3 | 1 | 0.050 s, degree 4 (success) | 0.016 s, degree 4 (success) | 0.008 s, degree 4 (success) | 0.012 s, degree 4 (success) |
| 1/3 | 3 | — s, degree — (no_certificate) | — s, degree — (no_certificate) | — s, degree — (no_certificate) | — s, degree — (no_certificate) |

## Reproduce

```sh
python3 tests/run_stage_division_benchmarking.py
```

The runner reads the fixed flags from `make -Bn` and independently toggles only the two requested flags.

Raw logs: `/var/folders/8d/y90177gs63xdbkbd34sv5wxr0000gn/T/stage_division_benchmarking_9i_4t_if`. Per-run data: `benchmarking_factor_stage_runs.json`.
