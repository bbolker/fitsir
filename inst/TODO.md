To do
==========

Lots and lots.

## Near future/low-hanging fruit

* finish up self-starting recipes
* implement fitting based on incidence rather than prevalence
* allow neg binomial rather than least-squares fit? log-Normal fit?
* plotting methods for data, fits?

## Optimization-related

* larger random/Latin hypercube sim; k-NN clustering, basins of attraction, etc.
* stochastic global opt? (Mullen JSS paper) 

## Longer-term

* add other data sets (*simple* univariate time series -- Harbin plague; etc.)
* think a lot about the interface: how to optimize efficiency, flexibility, breadth, transparency, simplicity ... ?  Want simple model objects (S3 class) that allow simple summaries, diagnostics, fits ... extend `mle2` class?
* should/must `mle2` class be extended via S4? or export `coef`?
* (`summary.pars` gets Roxygenized documentation confused)
* more general parameter transformation machinery -- from constrained to unconstrained scale to epidemiological parameters (R0, r, rho, etc.) and reverse
* flexibility in models?
* throw in old BIRS `erlangSEIR` stuff?
* gradient matching?
* state-space models: `pomp`, `JAGS`, etc. ?
* tie-in with Tiff Bogich's "Plom" stuff?
* tie-in with other epi R packages: `EpiModel`, `EpiDynamics` ...
