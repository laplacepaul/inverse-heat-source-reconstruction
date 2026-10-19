# Research documentation

This folder supplements the [project README](../README.md), the [nine notebooks](../notebooks/), and the [technical report](../report/). It is intended to make the numerical experiments independently understandable and their outputs inspectable.

| Document | What it explains |
|---|---|
| [METHODOLOGY.md](METHODOLOGY.md) | Model assumptions, Taylor/FTCS method, discrete inverse map, objective function and adjoint conventions |
| [EXPERIMENT_GUIDE.md](EXPERIMENT_GUIDE.md) | Aim, configuration, outputs and interpretation of each of nine notebooks |
| [RESULTS_SUMMARY.md](RESULTS_SUMMARY.md) | Audited numerical findings, comparisons and important qualifications |
| [DATA_DICTIONARY.md](DATA_DICTIONARY.md) | Every preserved data artifact, keys, array shapes and safe loading examples |
| [REPRODUCIBILITY.md](REPRODUCIBILITY.md) | Installation, notebook order, working-directory dependencies, verification and limitations |

The report is a technical manuscript, not a peer-reviewed paper. All experimental quantities refer to synthetic heat-equation data. A successful installation or the presence of saved outputs is not proof that every notebook was independently rerun in a clean environment.
