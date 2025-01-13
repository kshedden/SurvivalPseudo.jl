# SurvivalPseudo.jl

This package calculates fast [pseudo-observations](https://pubmed.ncbi.nlm.nih.gov/19654170/) for survival analysis in Julia.

We implement the [infinitessimal jackknife](https://cran.r-project.org/web/packages/survivalVignettes/vignettes/pseudo.html) and
an approach based on [martingale residuals](https://arxiv.org/abs/2109.02959).  These methods provide pseudo-observations for the
survival function.

We also include a function for calculating pseudo-observations for the [restricted mean life](https://link.springer.com/article/10.1007/s10985-004-4771-0).

See the examples folder for use.
