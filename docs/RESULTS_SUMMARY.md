# Verified numerical results and interpretation

The quantities below are derived from the saved experiment data or their companion JSON metrics. Units and error definitions matter: Notebook 05 reports **coefficient** errors, while the source reconstruction comparisons elsewhere commonly report **source-field** errors.

## A. Forward accuracy

The smooth analytical benchmark has an exact sinusoidal temperature solution. Notebook 01's joint refinement of spatial grid and stable timestep had recorded observed orders about 1.875 and 2.003. Notebook 09 independently recorded spatial orders **1.997927** and **1.999482**, and temporal orders **1.002749**, **1.001372**, **1.000685** (32, 64, 128, 256 steps respectively). The temporal verification compares FTCS against the exact **semi-discrete** solution; the spatial verification compares exact semi-discrete and exact continuous solutions.

## B. Three-coefficient reconstruction

Notebook 03 produces $A\in\mathbb{R}^{48	imes3}$ from 16 spatial sensors sampled at three times. It has full column rank and reported condition number **1.806**. Notebook 04's saved reconstructed source gives **0.5684% relative source-field error** and **0.6262% relative coefficient error** for the ridge solution. The output key `d_predicted` represents $A\widehat{	heta}$.

## C. Robustness: 150 realizations per setting

**Noise sweep:** Values are mean *coefficient* errors, in percent. The prescribed noise fraction is defined using the maximum clean observation magnitude; it is not necessarily the realized relative Euclidean noise norm.

| Prescribed noise fraction | Mean LS coefficient error | Mean ridge coefficient error |
|---:|---:|---:|
| 0 | 0.000% | 0.137% |
| 0.01 | 0.657% | 0.669% |
| 0.03 | 1.923% | 1.921% |
| 0.05 | 3.463% | 3.469% |
| 0.1 | 6.816% | 6.806% |

**Sensor sweep:** Values are mean *coefficient* errors over 150 draws for each sensor setting.

| Spatial sensors | Condition number | Mean LS coefficient error | Mean ridge coefficient error |
|---:|---:|---:|---:|
| 3 | 120.68 | 322.68% | 60.16% |
| 4 | 27.13 | 65.23% | 52.47% |
| 6 | 35.83 | 57.83% | 40.78% |
| 8 | 2.00 | 4.45% | 4.44% |
| 12 | 1.38 | 3.65% | 3.65% |
| 16 | 1.81 | 3.03% | 3.03% |

![Robustness experiments](../figures/phase5_robustness_experiments.png)

## D. Large and mismatched inverse systems

Notebook 06's 25-source-basis full system contains 48 observations and has a reported condition number near **1025**. The 18-observation version has rank at most 18 and a null-space dimension of at least seven (recorded rank 18). Regularization reduces errors in the tested noisy scenarios, but neither a unique ridge solution nor a small data residual establishes exact source identification.

Notebook 07 compares matched/representable source configurations against a source containing an omitted localized feature. Across grids $N=18,24,36,48$, source-field errors are about **4.1–4.2%** for representable noisy cases and about **8%** under noisy source-model mismatch. This gap is not removed by spatial refinement in the tested range.

![Discretization and model mismatch](../figures/phase7_discretization_model_mismatch.png)

## E. Adjoint verification

The saved Notebook 08 metrics record:

| Validation or reconstruction metric | Result |
|---|---:|
| Relative adjoint vs central finite-difference gradient discrepancy | 5.515757e-11 |
| Relative adjoint vs explicit matrix gradient discrepancy | 7.980403e-16 |
| Relative optimized vs direct coefficient discrepancy | 5.090569e-05 |
| Relative source-field error (25 coefficients) | 9.430% |
| Adjoint gradient evaluation (100 coefficients) | 0.0311 s |
| Explicit response matrix build (100 coefficients) | 0.8758 s |

The last two rows time different operations and do **not** establish an end-to-end speedup.

![Discrete adjoint validation](../figures/phase8_adjoint_results.png)

## F. Independent sensor-layout and end-to-end tests

Notebook 09 recorded 90 layouts per sensor-count group. Mean source-field error was **42.15%** for six sensors versus **20.96%** for sixteen sensors, with dependence on layout as well as number. The corrected benchmark recorded complete reconstruction medians of **0.126 s matrix vs 1.582 s adjoint** for 25 coefficients and **0.445 s matrix vs 2.810 s adjoint** for 100 coefficients. In these particular small linear problems, direct matrix-based reconstruction was faster.

![Independent verification and timings](../figures/task9_validation_figure.png)

## Boundaries of interpretation

- These are synthetic-data studies for a 2D linear heat model, not measurements from a physical experiment.
- The notebooks use different settings across experiments; differences in error percentages cannot automatically be assigned to one changed variable.
- A valid gradient is not equivalent to a uniquely identifiable source.
- Timing values depend on the execution environment; no comprehensive hardware-normalized benchmark is claimed.
- Consult the saved `.npz` and `.json` files and the corresponding notebook cells for detailed configurations.
