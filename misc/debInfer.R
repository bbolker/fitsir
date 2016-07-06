library(deSolve)
library(fitsir)

bombay2 <- setNames(bombay, c("tvec", "count"))

SIR.grad <- function(t, x, params) {
	g <- with(as.list(c(x,params)),{
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

library(deBInfer)
beta <- debinfer_par(name = "beta", var.type = "de", fixed = FALSE,
									value = 10, prior="unif", hypers=list(min = 0, max = 30),
									prop.var=3, samp.type="rw")

gamma <- debinfer_par(name = "gamma", var.type = "de", fixed = FALSE,
									value = 10, prior="unif", hypers=list(min = 0, max = 30),
									prop.var=3, samp.type="rw")


S <- debinfer_par(name = "S", var.type = "init", fixed = FALSE,
									value = 10000, prior="lnorm", hypers=list(meanlog = 10, sdlog = 1),
									prop.var=1000, samp.type="rw")

I <- debinfer_par(name = "I", var.type = "init", fixed = FALSE,
									value = 4, prior="unif", hypers=list(min = 0, max = 2*bombay2$count[1]),
									prop.var=0.1, samp.type="rw")

R <- debinfer_par(name = "R", var.type = "init", fixed = TRUE,
									value = 0)

mcmc.pars <- setup_debinfer(beta, gamma, I, S, R)

iter = 5000

mcmc_samples <- de_mcmc(N = iter, data=bombay2, de.model=SIR.grad,
												obs.model=SIR.logLik, all.params=mcmc.pars,
												Tmax = max(bombay2$tvec), data.times=bombay2$tvec, cnt=500,
												plot=TRUE, solver="ode")

plot(mcmc_samples)

burnin = 500
pairs(mcmc_samples, burnin = burnin, scatter=TRUE, trend=TRUE)

post_traj <- post_sim(mcmc_samples, n=500, times=1:31, burnin=burnin, output = 'all', prob = 0.95)

plot(post_traj, plot.type = "medianHDI", lty = c(2,1), lwd=3, col=c("red","grey20"),
		 panel.first=lines(bombay2, col="darkblue", lwd=2))

legend("topleft", legend=c("posterior median", "95% HDI", "true model"),
			 lty=c(2,1,1), lwd=c(3,2,2), col=c("red","grey20","darkblue"), bty='n')

