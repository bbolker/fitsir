library(deSolve)
library(fitsir)

bombay2 <- setNames(bombay, c("tvec", "count"))

SIR.grad <- function(t, x, params) {
	g <- with(as.list(c(x,params)),{
    
    gamma = r/R0.1
    beta = gamma + r
    
		dS = -beta*S*I/(S+I+R)
		dI = beta*S*I/(S+I+R)-gamma*I
		dR = gamma*I
		list(c(dI,dS,dR))
	})
}

dnorm2 <- function(x,mean,log=FALSE) {
	rmse <- sqrt(sum((x-mean))^2/(length(x)-1))
	return(dnorm(x,mean,sd=rmse,log=log))
}

SIR.logLik <- function(data, sim.data, samp){
	r <- sum(dnorm2(data$count, sim.data[,"I"], log = TRUE))
	return(r)
}

fitsir.deBInfer <- function(data,
                            iter = 10000,
                            burnin = 1000,
                            prob = 0.95,
                            start = NULL){
  tvec = data$tvec
  count = data$count
  
  is.null(start){
    start = list(R0 = 2.3, r = 0.6, S0 = 10000, I0 = 4)
  }
  
  with(as.list(start),{
    ## It's actually R0-1
    R0.1 <- debinfer_par(name = "R0.1", var.type = "de", fixed = FALSE,
                         value = R0-1, prior="unif", hypers=list(min = 0, max = 30),
                         prop.var=2, samp.type="rw")
    
    r <- debinfer_par(name = "r", var.type = "de", fixed = FALSE,
                      value = 0.6, prior="unif", hypers=list(min = 0, max = 30),
                      prop.var=0.3, samp.type="rw")
    
    S <- debinfer_par(name = "S", var.type = "init", fixed = FALSE,
                      value = S0, prior="unif", hypers=list(min = 100, max = 100000),
                      prop.var=10000, samp.type="rw")
    
    I <- debinfer_par(name = "I", var.type = "init", fixed = FALSE,
                      value = I0, prior="unif", hypers=list(min = 0, max = 2*bombay2$count[1]),
                      prop.var=2, samp.type="rw")
    
    R <- debinfer_par(name = "R", var.type = "init", fixed = TRUE,
                      value = 0)
    
    mcmc.pars <- setup_debinfer(R0.1, r, I, S, R)
    
    iter = 10000
    
    mcmc_samples <- de_mcmc(N = iter, data=bombay2, de.model=SIR.grad,
                            obs.model=SIR.logLik, all.params=mcmc.pars,
                            Tmax = max(bombay2$tvec), data.times=bombay2$tvec, cnt=500,
                            plot=TRUE, solver="ode")
    
  })
  
}

library(deBInfer)

plot(mcmc_samples)

burnin = 1000
pairs(mcmc_samples, burnin = burnin, scatter=TRUE, trend=TRUE)

post_traj <- post_sim(mcmc_samples, n=500, times=1:31, burnin=burnin, output = 'all', prob = 0.95)

plot(post_traj, plot.type = "medianHDI", lty = c(1,2), lwd=3, col=c("red","grey20"),
     log = "y",
		 panel.first=lines(bombay2, col="darkblue", lwd=2))

plot(bombay2, log = "y")
lines(post_traj[["median"]]$I)
matlines(post_traj[["HDI"]]$I, col = 2, lty = c(2,2))
