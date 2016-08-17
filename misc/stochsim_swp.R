if (file.exists("~/Rlibs")) .libPaths("~/Rlibs")

library(fitsir)
library(plyr)
library(devtools)
source("stochsim_funs.R")
findSens <- fitsir:::findSens
fitsir.optim <- fitsir:::fitsir.optim
SIR.logLik <- fitsir:::SIR.logLik
## beta_range <- c(0.16,0.75)
## sample over R0 and post-process to beta=R0/gamma, to
##   avoid R0<1 cases from sampling beta and gamma independently ...
R0_range <- c(1.1,8)
gamma_range <- c(0.07,0.3)
N_range <- c(500,3e4)
I0_range <- c(1,20)

nsim <- 500
lhs_df <- cbind(
    ## beta=seq(beta_range[1],beta_range[2],length.out=nsim),
    R0=seq(R0_range[1],R0_range[2],length.out=nsim),
    gamma=seq(gamma_range[1],gamma_range[2],length.out=nsim),
    N=seq(N_range[1],N_range[2],length.out=nsim),
    I0=seq(I0_range[1],I0_range[2],length.out=nsim))
set.seed(101)
for (i in 2:ncol(lhs_df)) {
    lhs_df[,i] <- sample(lhs_df[,i])
}
fn <- "stochsim_swp.rda"

g <- SIR.logLik()
truepar <- as.data.frame(matrix(NA, ncol = 4, nrow = nsim))
names(truepar) <- c("log.beta", "log.gamma", "log.N", "logit.i")

tmpfun <- function(i, fitfun = fitsir, ...) {
    p <- lhs_df[i,]
    p2 <- c(beta=unname(p["R0"]*p["gamma"]),
            p[2:3],
            i0=unname(p["I0"]/p["N"]))
    true.pars <- with(as.list(p2),
               c(log.beta=log(beta),log.gamma=log(gamma),
                 log.N=log(N),logit.i=qlogis(i0)))
    d <- simfun(pars=p2,tmax=100,dt=1,rpars=list(size=3),seed=101)
    true.nll <- c("true_nll" = g(true.pars, d$count))
    truefit <- fitfun(d, start = true.pars, ...)
    if(class(truefit) != "numeric"){
        truefit.pars <- coef(truefit)
    }else{
        truefit.pars <- truefit
    }
    truefit.nll <- c("truefit_nll" = g(truefit.pars, d$count))
    start.pars <- unlist(startfun(data = d, auto = TRUE))
    start.nll <- c("start_nll" = g(start.pars, d$count))
    fit <- fitfun(d, start = start.pars, ...)
    if(class(truefit) != "numeric"){
        fit.pars <- coef(fit)
    }else{
        fit.pars <- fit
    }
    fit.nll <- c("fit_nll" = g(fit.pars, d$count))
    names(true.pars) <- paste0("true_", names(true.pars))
    names(truefit.pars) <- paste0("truefit_", names(truefit.pars))
    names(start.pars) <- paste0("start_", names(start.pars))
    names(fit.pars) <- paste0("fit_", names(start.pars))
    return(c(true.pars, true.nll, truefit.pars, truefit.nll, start.pars, start.nll, fit.pars, fit.nll))
}

system.time(exfit <- tmpfun(2, fitfun = fitsir)) ##13.36
system.time(tmpfun(2, fitfun = fitsir.optim)) ##4.46
system.time(tmpfun(2, fitfun = fitsir.optim, nll = TRUE)) ##5.22

res.fitsir <- res.fitsir.optim <- res.fitsir.optim.nll <- as.data.frame(matrix(NA, ncol = length(exfit), nrow = nsim))

names(res.fitsir) <- names(res.fitsir.optim) <- names(res.fitsir.optim.nll) <- names(exfit)

for(i in 1:nsim){
    cat(i, "\n")
    f <- try(tmpfun(i))
    cat("f\n")
    f.optim <- try(tmpfun(i, fitfun = fitsir.optim))
    cat("f.optim\n")
    f.optim2 <- try(tmpfun(i, fitfun = fitsir.optim, nll = TRUE))
    cat("f.optim2\n")
    if (!is(f,"try-error")) {
        res.fitsir[i,] <- f
    }
    if (!is(f.optim,"try-error")) {
        res.fitsir.optim[i,] <- f.optim
    }
    if (!is(f.optim2,"try-error")) {
        res.fitsir.optim.nll[i,] <- f.optim2
    }
    save("res.fitsir", "res.fitsir.optim", "res.fitsir.optim.nll",file=fn)
}