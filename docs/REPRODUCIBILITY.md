# Reproducing the numerical experiments

## Status and scope

These notebooks were previously executed in Google Colab, and saved output files are provided in `data/`. The repository has **not** been certified as a clean-kernel, clean-environment reproduction of all nine notebooks. Python package versions and execution hardware were not fully pinned at the time of the original runs.

## Setup

From the project root, use Python 3 and a virtual environment:

```bash
python -m venv .venv
source .venv/bin/activate  # on Windows: .venv\Scripts\activate
python -m pip install -r requirements.txt
python -m pip install jupyterlab
jupyter lab
```

Alternatively, open the numbered notebooks in Google Colab. The notebook working directory may differ from the repository root in Colab. Always check file input/output paths before execution.

## Working-directory dependencies

- Start with notebooks 01–02 to verify the forward solver and spatial sources.
- Notebook 03 **writes** `phase3_observation_data.npz` into the active notebook working directory.
- Notebooks 04 and 05 **read** `phase3_observation_data.npz` by filename. Place a copy in their current working directory before execution, or change the path to `../data/phase3_observation_data.npz` when running from `notebooks/`.
- Notebooks 06–09 principally assemble their experiment models internally and write their own named result files. They should still be inspected for Colab-specific paths, dataset saving, and large computations.
- Saved outputs are available in `data/` for inspection without rerunning experiments. Figures are in `figures/`.

**Important:** Saved result files are not proof that rerunning a notebook with a different random seed or software version will produce identical numbers. Fix random seeds, record software versions, and run with a controlled execution environment to make numerical comparisons meaningful.

## Expected validation checkpoints

| Check | Expected behaviour |
|---|---|
| Zero source, zero initial state | Numerically zero temperature |
| Homogeneous Dirichlet boundaries | Fixed zero boundary values |
| FTCS timestep | $r_x+r_y\le 1/2$ for every forward run |
| Small response matrix | Shape `(48, 3)` with full column rank |
| Gradient validation | Small finite-difference and matrix-gradient discrepancies |
| Notebook 09 temporal steps | `[32, 64, 128, 256]` only |
| Runtime comparisons | Do not compare one gradient evaluation against a complete direct reconstruction |

The original 16-step temporal configuration was removed from Notebook 09's corrected analysis because it violated the standard FTCS stability restriction.

## Check packaged artifacts without executing notebooks

Run from the repository root:

```python
from pathlib import Path
import json
import numpy as np

root = Path("data")
with np.load(root / "phase3_observation_data.npz", allow_pickle=False) as f:
    assert f["A"].shape == (48, 3)

with np.load(root / "task9_validation_data.npz", allow_pickle=False) as f:
    print("Saved keys:", f.files)

m = json.loads((root / "phase8_metrics.json").read_text())
print("Adjoint/matrix gradient check:",
      m["adjoint_vs_matrix_gradient_relative_error"])
```

## What a full reproducibility certification would require

1. Record tested Python and dependency versions and platform details.
2. Execute each notebook from a clean kernel in its documented folder or a consistent root directory.
3. Remove implicit dependence on Colab's temporary filesystem, making paths relative to the repository.
4. Record seeds for randomized tests and keep output metrics alongside generation code.
5. Compare regenerated figures, errors, ranks, convergence slopes and benchmark protocols against this release; expect hardware-dependent differences in timing.

See [EXPERIMENT_GUIDE.md](EXPERIMENT_GUIDE.md) for research questions and [DATA_DICTIONARY.md](DATA_DICTIONARY.md) for artifact contents.
