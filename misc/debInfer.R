library(deSolve)
library(deBInfer)
library(fitsir)
source("stochsim_funs.R")

data <- simfun(pars = c(beta = 2, gamma = 1, N = 10000, i0 = 0.01),
    rpars = list(sd = 50), seed = 101)

SIR.grad <- function(t, x, params) {
	g <- with(as.list(c(x,params)),{
        
	    gamma <- r/(R0 - 1)
	    beta <- gamma + r
	    
		dS <- -beta*S*I/(S+I+R)
		dI <- beta*S*I/(S+I+R)-gamma*I
		dR <- gamma*I
		list(c(dI,dS,dR))
	})
}

dnorm2 <- fitsir:::dnorm2

SIR.logLik <- function(data, sim.data, samp){
	r <- sum(dnorm2(data$count, sim.data[,"I"], log = TRUE))
	return(r)
}

fitsir.deBInfer <- function(data,
                            R0.args = list(),
                            r.args = list(),
                            S.args = list(),
                            I.args = list(),
                            iter = 10000,
                            plot = FALSE,
                            cnt = 500, ...){
    tvec <- data$tvec
    count <- data$count
    
    ## FIXME: Enforce DRY coding
    ## Can't think of an easy way...
    default <- list(fixed = FALSE, samp.type = "rw", prior = "unif")
    R0.default <- append(list(value = 2, hypers = list(min = 1, max = 10), prop.var = 2), default)
    r.default <- append(list(value = 1, hypers = list(min = 0.01, max = 4), prop.var = 0.1), default)
    S.default <- append(list(value = 2000, hypers = list(min = 1e3, max = 1e6), prop.var = 1000), default)
    I.default <- append(list(value = 10, hypers = list(min = 1, max = 20), prop.var = 2), default)
    
    defaultlist <- list(R0.default, r.default, S.default, I.default)
    
    argslist <- list(R0.args, r.args, S.args, I.args)
    
    for(i in 1:4){
        def <- defaultlist[[i]]
        new <- argslist[[i]]
        
        m <- match(names(new), names(def))
        
        defaultlist[[i]] <- replace(def, m, new)
    }
    
    ## FIXME: Enforce DRY coding here as well...
    R0 <- do.call("debinfer_par", args = append(defaultlist[[1]], list(name = "R0", var.type = "de")))
    r <- do.call("debinfer_par", args = append(defaultlist[[2]], list(name = "r", var.type = "de")))
    S <- do.call("debinfer_par", args = append(defaultlist[[3]], list(name = "S", var.type = "init")))
    I <- do.call("debinfer_par", args = append(defaultlist[[4]], list(name = "I", var.type = "init")))
    R <- debinfer_par(name = "R", var.type = "init", fixed = TRUE, value = 0)
    
    mcmc.pars <- setup_debinfer(R0, r, I, S, R)
        
    mcmc_samples <- de_mcmc(N = iter, data=data, de.model=SIR.grad,
        obs.model=SIR.logLik, all.params=mcmc.pars,
        Tmax = max(tvec), data.times=tvec, cnt=cnt,
        plot=plot, solver="ode", ...)
}

## sample run demonstrating how to use the arguments
f <- fitsir.deBInfer(data,
    R0.args = list(value = 3, prior = "lnorm", hypers = list(meanlog = 2, sdlog = 1)),
    I.args = list(value = 50, hypers = list(min =30, max = 100)),
    iter = 2000,
    plot = TRUE)

plot(f)

burnin <- 500
pairs(f, burnin = burnin, scatter=TRUE, trend=TRUE)

post_traj <- post_sim(f, n=500, times=data$tvec, burnin=burnin, output = 'all', prob = 0.95)

plot(post_traj, plot.type = "medianHDI", lty = c(1,2), lwd=3, col=c("red","grey20"),
		 panel.first=lines(data, col="darkblue", lwd=2))

plot(data, log = "y")
lines(data$tvec, post_traj[["median"]]$I)
matlines(data$tvec, post_traj[["HDI"]]$I, col = 2, lty = c(2,2))
