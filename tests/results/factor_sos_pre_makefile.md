# FACTOR_SOS_PRE comparison with current Makefile settings

This repeat changes only `FACTOR_SOS_PRE`. Both runs use the current Makefile configuration: `ENABLE_BINARY_SEARCH_AVKL`, `ENABLE_DEBUGGING`, `EPS_FACTOR=17/20`, `N_GUESS_AVKL=70`, and `N_GUESS_LS=600`. Weifeng and all other optional algorithm flags remain disabled. The same benchmark source, cases, validation, output redirection, and measurement procedure are used in both runs. The `FACTOR_SOS_BENCHMARK` harness guard is defined in both. No Makefile or installed .mla changes were made.

The regression test's constants now act as defaults, so command-line values are honored rather than overwritten by its former hardcoded N guesses of 2/2.

## Method

Maple 2022; sequential disabled and enabled runs. Source loading is excluded. Each operation is warmed up, batches are calibrated to at least 50 ms or at most 4096 calls, and seven subsequent wall-time batches are measured. Values below are median milliseconds per call. The full regression suite includes validation; individual construction rows exclude validation. Debug output remains enabled in both runs and is redirected to log files. Timing differences, especially on small cases and early returns, have substantial variability and should not be interpreted as statistically established effects.

Certificate degree vectors describe the expanded coefficients `[sigma_0, sigma_1, ...]`. A zero coefficient is shown as `zero`. Maximum certificate degree below means maximum coefficient degree, not degree of a square summand. For the real pipeline cases, it also equals the maximum degree of weighted terms `sigma_i*basis_i`. For the mocked controlled stages, coefficient degree is 6 and weighted degree is 8, unchanged by the flag.

## Results

Negative changes mean faster with factoring enabled.

| Case | Off (ms) | On (ms) | Difference (ms) | Change (%) | Coefficient degrees off → on |
|---|---:|---:|---:|---:|---|
| regression_suite | 360.000 | 308.000 | -52.000 | -14.4% | — |
| isSOS_positive | 1.344 | 0.375 | -0.969 | -72.1% | — |
| isSOS_sign_change | 0.003 | 0.500 | 0.497 | 18515.0% | — |
| interval_positive | 5.625 | 0.938 | -4.688 | -83.3% | — |
| interval_sign_change | 1.469 | 2.000 | 0.531 | 36.2% | — |
| minimum | 4.625 | 3.813 | -0.813 | -17.6% | — |
| epsilon | 24.750 | 25.500 | 0.750 | 3.0% | — |
| lower_bound_constant | 2.563 | 2.906 | 0.344 | 13.4% | — |
| lower_bound | 38.000 | 40.000 | 2.000 | 5.3% | — |
| controlled_stages | 12.750 | 15.250 | 2.500 | 19.6% | [6, 6, 6] → [6, 6, 6] |
| certificate_factored_input | 84.000 | 28.500 | -55.500 | -66.1% | [18, 16] → [4, 2] |
| certificate_linear | 25.250 | 29.000 | 3.750 | 14.9% | [2, 0] → [2, 0] |
| certificate_global_SOS | 4.625 | 4.500 | -0.125 | -2.7% | [4, zero] → [4, zero] |
| certificate_minus_one | 4.813 | 7.250 | 2.438 | 50.6% | [6, 6] → [6, 6] |

The main factored-input case, f=(x^2+1)*(x+2), basis=[1-x^2], took **84.00 → 28.50 ms**, a **55.50 ms reduction (66.1%)**, and its maximum certificate degree fell **18 → 4**. Its observed batch ranges were 80–212 ms without factoring and 27.5–49 ms with factoring. The full regression suite took **360 → 308 ms**, a **52 ms reduction (14.4%)**. These are observed timings, not a guarantee of speedup on other inputs.

The linear case f=x+1, basis=[-x^2], retained degree 2. The globally SOS input f=(x^2+1)*(x^2+2), basis=[-x^2], retained degree 4; f=(x^2+1)*(x+1), basis=[-1], retained degree 6. Both early-return cases execute before the heuristic, so their timing differences reflect measurement variability rather than factoring work. The controlled-stage example deliberately has identical final coefficients and uses mocked construction steps; it checks factor bookkeeping rather than end-to-end performance.

All regression checks passed in both runs, including exact certificate reconstruction and nonnegativity. Full per-case timing ranges and calibrated call counts are in `factor_sos_pre_makefile_timings.txt`.

## Exact commands

```sh
maple -q -DENABLE_BINARY_SEARCH_AVKL -DENABLE_DEBUGGING -DEPS_FACTOR=17/20 -DN_GUESS_AVKL=70 -DN_GUESS_LS=600 tests/benchmark_factor_sos_pre.mpl > /tmp/factor_sos_current_off.log
maple -q -DENABLE_BINARY_SEARCH_AVKL -DENABLE_DEBUGGING -DEPS_FACTOR=17/20 -DN_GUESS_AVKL=70 -DN_GUESS_LS=600 -DFACTOR_SOS_PRE tests/benchmark_factor_sos_pre.mpl > /tmp/factor_sos_current_on.log
```

