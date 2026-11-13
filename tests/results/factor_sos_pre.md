# SOS factor preprocessing benchmark

Measured with Maple 2022 using the cases in `tests/factor_sos_pre.mpl`.

## Method

`tests/benchmark_factor_sos_pre.mpl` loads the source directly and validates the regression suite before timing. Each operation is warmed up, batch size doubles until it takes at least 50 ms (or reaches 4096 calls), and seven subsequent batches are measured with `time[real]()`. Tables show the median wall time per call in milliseconds. Configurations run sequentially. Source loading is excluded; the `regression_suite` row includes all validation, while individual certificate rows time construction only. EPS_FACTOR=17/20, N_GUESS_AVKL=2, and N_GUESS_LS=2 match the tests. No other optional algorithm flags were enabled.

Degree vectors list the exact expanded polynomial degrees of `[sigma_0, sigma_1, ...]` in `f = sigma_0 + sum(sigma_i*basis_i)`. `zero` denotes a zero coefficient (Maple reports degree -infinity). These are polynomial coefficient degrees, not degrees of individual square summands. For the four real pipeline cases, the maximum weighted degree `max(degree(sigma_0), degree(sigma_i*basis_i))` equals the maximum coefficient degree shown. For the constructed-stage case, maximum coefficient degree is 6 and maximum weighted degree is 8 in both configurations.

Time change is `(enabled/disabled - 1)*100`; negative means faster. These small cases have substantial timing variability. Raw batch ranges are saved in `factor_sos_pre_timings.txt`; small differences should not be treated as statistically established improvements.

## WEIFENG_OPTIMIZATION disabled

| Case | Off (ms) | On (ms) | Change (ms) | Change (%) | Coefficient degrees off → on |
|---|---:|---:|---:|---:|---|
| regression_suite | 282.000 | 213.000 | -69.000 | -24.5% | — |
| isSOS_positive | 1.359 | 0.430 | -0.930 | -68.4% | — |
| isSOS_sign_change | 0.005 | 0.492 | 0.488 | 10509.8% | — |
| interval_positive | 1.438 | 0.285 | -1.152 | -80.2% | — |
| interval_sign_change | 0.523 | 1.313 | 0.789 | 150.7% | — |
| minimum | 1.688 | 0.891 | -0.797 | -47.2% | — |
| epsilon | 21.750 | 18.000 | -3.750 | -17.2% | — |
| lower_bound_constant | 2.188 | 2.188 | 0.000 | 0.0% | — |
| lower_bound | 29.500 | 32.500 | 3.000 | 10.2% | — |
| controlled_stages | 8.250 | 10.500 | 2.250 | 27.3% | [6, 6, 6] → [6, 6, 6] |
| certificate_factored_input | 65.000 | 18.250 | -46.750 | -71.9% | [10, 8] → [4, 2] |
| certificate_linear | 14.250 | 14.500 | 0.250 | 1.8% | [2, 0] → [2, 0] |
| certificate_global_SOS | 4.469 | 4.438 | -0.031 | -0.7% | [4, zero] → [4, zero] |
| certificate_minus_one | 4.813 | 4.688 | -0.125 | -2.6% | [6, 6] → [6, 6] |

## WEIFENG_OPTIMIZATION enabled

| Case | Off (ms) | On (ms) | Change (ms) | Change (%) | Coefficient degrees off → on |
|---|---:|---:|---:|---:|---|
| regression_suite | 311.000 | 318.000 | 7.000 | 2.3% | — |
| isSOS_positive | 1.359 | 0.414 | -0.945 | -69.5% | — |
| isSOS_sign_change | 0.004 | 0.547 | 0.542 | 12343.1% | — |
| interval_positive | 1.438 | 0.285 | -1.152 | -80.2% | — |
| interval_sign_change | 0.828 | 1.016 | 0.188 | 22.6% | — |
| minimum | 1.625 | 0.891 | -0.734 | -45.2% | — |
| epsilon | 20.000 | 22.250 | 2.250 | 11.3% | — |
| lower_bound_constant | 1.781 | 2.188 | 0.406 | 22.8% | — |
| lower_bound | 24.500 | 33.500 | 9.000 | 36.7% | — |
| controlled_stages | 8.000 | 10.125 | 2.125 | 26.6% | [6, 6, 6] → [6, 6, 6] |
| certificate_factored_input | 51.000 | 59.000 | 8.000 | 15.7% | [6, 4] → [12, 10] |
| certificate_linear | 55.000 | 55.000 | 0.000 | 0.0% | [10, 8] → [10, 8] |
| certificate_global_SOS | 4.563 | 4.438 | -0.125 | -2.7% | [4, zero] → [4, zero] |
| certificate_minus_one | 4.563 | 4.563 | 0.000 | 0.0% | [6, 6] → [6, 6] |

## Inputs and interpretation

- `certificate_factored_input`: f=(x^2+1)*(x+2), basis=[1-x^2].
- `certificate_linear`: f=x+1, basis=[-x^2].
- `certificate_global_SOS`: f=(x^2+1)*(x^2+2), basis=[-x^2]. This takes the global nonnegativity early return before preprocessing.
- `certificate_minus_one`: f=(x^2+1)*(x+1), basis=[-1]. This also returns before preprocessing.
- `controlled_stages` uses mocked construction steps to validate cumulative factor restoration. Its degree equality is intentional; its runtime is not representative of the full algorithm. Its input differs between Weifeng modes as specified in the tests.

Without Weifeng, the factored-input case improved from 65.00 to 18.25 ms (71.9% faster), with maximum certificate degree reduced from 10 to 4. With Weifeng, the same case worsened from 51.00 to 59.00 ms (15.7% slower), with maximum degree increasing from 6 to 12. Reducing the input changes the constants and exponents selected by the downstream heuristic, so it does not guarantee lower final certificate degree. The other real pipeline cases had unchanged degrees. All reported configurations passed the regression suite, including exact certificate reconstruction and coefficient nonnegativity.

The sign-changing isSOS example exits immediately on a negative value at x=0 when preprocessing is disabled; factoring first adds overhead. Positive root-free and interval-positive examples show the clearest predicate-level improvements. The full suite was 24.5% faster without Weifeng and 2.3% slower with Weifeng in these measurements.

## Reproduce

```sh
maple -q tests/benchmark_factor_sos_pre.mpl
maple -q -DFACTOR_SOS_PRE tests/benchmark_factor_sos_pre.mpl
maple -q -DWEIFENG_OPTIMIZATION tests/benchmark_factor_sos_pre.mpl
maple -q -DWEIFENG_OPTIMIZATION -DFACTOR_SOS_PRE tests/benchmark_factor_sos_pre.mpl
```

