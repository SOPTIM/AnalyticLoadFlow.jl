# AnalyticLoadFlow workshop notebooks

The notebooks run on Google Colab, no local installation required. The first
cell installs AnalyticLoadFlow from GitHub (branch `main`) into a fresh
temporary environment (takes a few minutes); a commented line in the same
cell switches to the latest registered release or to a private checkout.

| Notebook | What it covers | Open |
|---|---|---|
| [workshop_apslf.ipynb](workshop_apslf.ipynb) | Part 1: the two-bus and four-bus hand calculations of theory Section 7 checked digit by digit against the solver, the old recursion without the conjugation for comparison, the PV kernels on the same network. Part 2: the 9-bus case, series coefficients and Padé, the stability indicator, the germ variants, PV buses and reactive limits, the sparse path | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/SOPTIM/AnalyticLoadFlow.jl/blob/main/notebooks/workshop_apslf.ipynb) |
| [workshop_pst.ipynb](workshop_pst.ipynb) | Transformer branch model, why a PST breaks the flat germ, the two embeddings of theory Section 6.5 worked by hand on a 4-bus network, PST angle and ratio sweeps on the 9-bus ring, a regulated PST via the outer loop | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/SOPTIM/AnalyticLoadFlow.jl/blob/main/notebooks/workshop_pst.ipynb) |
| [workshop_large_network.ipynb](workshop_large_network.ipynb) | PEGASE 2869 from a MATPOWER file: import with the MATPOWER conventions, sparse solve and timing, deviation vs no-load embedding on a large network, reactive limits, branch flows | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/SOPTIM/AnalyticLoadFlow.jl/blob/main/notebooks/workshop_large_network.ipynb) |

**Do not edit the `.ipynb` files directly**: they are generated. Edit the
Literate.jl source in [`docs/lit/`](../docs/lit/) and regenerate with
`julia --project=docs docs/generate_notebooks.jl`.

`julia --project=. docs/test_notebooks.jl` runs every source as a script in a
fresh Julia process and fails on an error or a warning; the `@assert ... #src`
lines in the sources (dropped from the notebooks) check the results against
the text. CI runs it on every push.
