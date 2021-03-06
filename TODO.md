# Overall to-do list



## Package stuff

- self-starting
    - prevalence is OK, but a little hacky (not enabled by default)
	- incidence is not so good
- work on neg binomial fit ... quasi-Poisson
- plotting methods for data, fits?

## Optimization-related

- larger random/Latin hypercube sim; k-NN clustering, basins of attraction, etc.
- stochastic global opt? (Mullen JSS paper)

## Other data sets

- 1918 flu

## Longer-term/lower priority

- priors
- missing data


- think a lot about the interface: how to optimize efficiency, flexibility, breadth, transparency, simplicity ... ?  Want simple model objects (S3 class) that allow simple summaries, diagnostics, fits ... extend `mle2` class?
- should/must `mle2` class be extended via S4? or export `coef`?
- (`summary.pars` gets Roxygenized documentation confused)
- more general parameter transformation machinery -- from constrained to unconstrained scale to epidemiological parameters (R0, r, rho, etc.) and reverse
- flexibility in models?
- throw in old BIRS `erlangSEIR` stuff?
- gradient matching?
- state-space models: `pomp`, `JAGS`, etc. ?
- tie-in with Tiff Bogich's "Public Library of Models" ([plom.io](https://github.com/tiffbogich/PLoM.io)) (related to Ballesteros's [epi-models](https://github.com/sballesteros/epi-models))?
- tie-in with other epi R packages: `EpiModel`, `EpiDynamics`, `pomp` ...

## exploration to-do list

- Figure 1: prevalence?? or preferably re-do with incidence if not too hard
   - (a) cumulative distributions of attained log-likelihood for N-M and BFGS, showing three levels (pick an example with three)
   - (b) points and trajectory for each plateau (i.e., flat, jump-up plus decay, good fit)
   - (c) ??? basins of attraction for difference plateaus: i.e. plot the starting points from the LHS coloured accord to which plateau they end up in.  Summaries???
- check Hessians (analytical + finite-diff of grad, which you can get via numDeriv)

- do this for real epidemic curves
- re-run resList.grad (BFGS) with negative log-likelihood outcomes
- big LHS experiment:
    - sensitivity equations vs. no explicit gradient
	- different optimizers (via optimx?)
	    - no point in running derivative-free methods with sens equations
		- in optimx, can't run no-gradient with derivative-ful methods
		- ?check latest optim (optimrx from r-forge?)
		   install.packages("optimrx",repos="http://r-forge.r-project.org")
    - NLL vs least squares?
	- rk4 vs lsoda
- compare:
    - simulated data
	- real data
	    - Bombay
		- Harbin
		- Phila 1918 flu
		- Gani & Leach smallpox?
		- ???
i.e., is it just that the Bombay data are weird?

- visualization

- try with 2 parameters (fix gamma, N, show likelihood
   surface for beta, I0); 3 parameters?
   
- reparameterization for most independent results?

- understand why NLL vs least-squares should matter.
    - just numerical sensitivity/"noise"?
	- should get same *best fit* value (extreme points of surface should agree ...)
