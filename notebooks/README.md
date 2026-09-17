# AnalyticLoadFlow workshop notebook

The notebook runs on Google Colab, no local installation required. The first
cell installs AnalyticLoadFlow from GitHub (branch `main`) into a fresh
temporary environment (takes a few minutes); a commented line in the same
cell switches to the latest registered release or to a private checkout.

| Notebook | What it covers | Open |
|---|---|---|
| [workshop_apslf.ipynb](workshop_apslf.ipynb) | Part 1: the two-bus and four-bus hand calculations of theory Section 7 checked digit by digit against the solver, the old recursion without the conjugation for comparison, the PV kernels on the same network. Part 2: the 9-bus case, series coefficients and Padé, the stability indicator, the germ variants, PV buses and reactive limits, the sparse path | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/SOPTIM/AnalyticLoadFlow.jl/blob/main/notebooks/workshop_apslf.ipynb) |

**Do not edit the `.ipynb` file directly**: it is generated. Edit the
Literate.jl source in [`docs/lit/`](../docs/lit/) and regenerate with
`julia --project=docs docs/generate_notebooks.jl`.
