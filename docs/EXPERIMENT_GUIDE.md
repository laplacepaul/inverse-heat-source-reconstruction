# Experiment-by-experiment guide

The experiments form a sequence, but some later notebooks construct their own numerical models. Notebook 03 is the key dependency for Notebooks 04 and 05. Open each notebook for exact experiment parameters; numerical parameters and noise draws can differ between experiments.

| Notebook | Research question | Key output or check | Saved artifacts |
|---|---|---|---|
| [01 — Forward solver](../notebooks/01_forward_heat_equation.ipynb) | Does FTCS solve the heat equation accurately? | Exact sinusoidal reference; 12/24/48 joint refinement | Results printed in notebook |
| [02 — Spatial sources](../notebooks/02_general_spatial_heat_sources.ipynb) | Can the solver handle Gaussian heat sources and superposition? | Gaussian fields; source linearity check | Results printed in notebook |
| [03 — Observation operator](../notebooks/03_observation_operator.ipynb) | Can responses to individual bases form a linear observation map? | 48-by-3 matrix, rank and condition; 16 sensors × 3 times | `phase3_observation_data.npz` |
| [04 — Inverse problem](../notebooks/04_inverse_problem.ipynb) | Can LS and ridge recover three coefficients from noisy temperatures? | Coefficient/field reconstruction; BFGS consistency | `phase4_reconstruction_results.npz` |
| [05 — Robustness](../notebooks/05_robustness.ipynb) | How do noise, sensors and penalty magnitude affect coefficients? | 150 Monte Carlo draws per noise/sensor setting | `phase5_experiment_results.npz`, Figure 1 |
| [06 — Larger inverse](../notebooks/06_large_inverse.ipynb) | How does source dimension affect identifiability? | 25 bases; 48 and 18 observations; regularization | `phase6_large_inverse_data.npz` |
| [07 — Mismatch](../notebooks/07_discretization_model_mismatch.ipynb) | Can grid refinement remove source-model mismatch? | N=18/24/36/48; representable vs unrepresented source features | `phase7_results.npz`, companion JSON, Figure 2 |
| [08 — Discrete adjoint](../notebooks/08_discrete_adjoint.ipynb) | Can we evaluate the parameter gradient without forming A? | Central-difference and matrix gradient checks; L-BFGS-B | `phase8_results.npz`, companion JSON, Figure 3 |
| [09 — Independent validation](../notebooks/09_validation.ipynb) | Do spatial/time orders and performance results survive independent tests? | Separate convergence, 90 randomized layouts per group, workflow runtime | `task9_validation_data.npz`, companion JSON, Figure 4 |

## Reading order

For numerical PDEs, begin at Notebooks 01–02 and Chapter 4 of the report. For the inverse problem, read Notebooks 03–07 together with Chapter 5. For adjoints and computational comparisons, use Notebooks 08–09 and Chapter 6–7.

## Interpreting the outputs

- Task 05 error scores are **coefficient errors**; do not substitute them for field errors.
- A matrix can have a finite ratio of positive singular values and still be non-identifiable when it has a nontrivial null space.
- Task 07 uses a source feature absent from the inversion dictionary; refining the grid cannot make that feature representable by the fixed dictionary.
- Task 08 compares **one gradient evaluation** with **one response-matrix build**. Task 09 compares **end-to-end reconstruction workflows**; the latter is the appropriate comparison for total timing.
- Task 09's corrected time-step sequence is 32, 64, 128, 256. An earlier 16-step experiment was unstable and is excluded.
