To do
==========

Lots and lots.

* think a lot about the interface: how to optimize efficiency, flexibility, breadth, transparency, simplicity ... ?  Want simple model objects (S3 class) that allow simple summaries, diagnostics, fits ... extend `mle2` class?
* should/must `mle2` class be extended via S4?
* (`summary.pars` gets Roxygenized documentation confused)
* more general parameter transformation machinery -- from constrained to unconstrained scale to epidemiological parameters (R0, r, etc.) and reverse
* plotting methods for data, fits?
* allow neg binomial rather than least-squares fit?
* implement fitting based on incidence rather than prevalence?
* flexibility in models?
* throw in old BIRS `erlangSEIR` stuff?
* gradient matching?
* state-space models: `pomp`, `JAGS`, etc. ?
* tie-in with Tiff Bogich's "Plom" stuff?
* tie-in with other epi R packages: `EpiModel`, `EpiDynamics` ...
* Bombay data fitting!
    * self-starting recipes
    * larger random/Latin hypercube sim; k-NN clustering, basins of attraction, etc.
    * stochastic global opt? (Mullen JSS paper) 
    * save vignette PDF; save
