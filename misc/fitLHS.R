library("deSolve")
library("fitsir")
library("devtools")
## this is for collywobbles to get away with the permission problem...
source("fitSIR_funs.R")
##
## load_all("..")
bombay2 <- setNames(bombay, c("tvec", "count"))

nsim <- 500

log10R0_range <- c(-3,1) ## log(R0-1)
log10gamma_range <- c(-1,2)
log10N_range <- c(2,8)
log10I0_range <- c(0,2)

lhs_df <- cbind(
    ## beta=seq(beta_range[1],beta_range[2],length.out=nsim),
    R0=1+10^(seq(log10R0_range[1],log10R0_range[2],length.out=nsim)),
    gamma=10^(seq(log10gamma_range[1],log10gamma_range[2],length.out=nsim)),
    N=10^(seq(log10N_range[1],log10N_range[2],length.out=nsim)),
    I0=10^(seq(log10I0_range[1],log10I0_range[2],length.out=nsim)))
set.seed(101)
for (i in 1:ncol(lhs_df)) {
    lhs_df[,i] <- sample(lhs_df[,i])
}

tmpf <- function(pars.vec){
    with(as.list(c(pars.vec)),{
        tpars <- c(
            log.beta = log(R0 * gamma),
            log.gamma = log(gamma),
            log.N = log(N),
            logit.i = qlogis(I0/N)
            
        )
        return(tpars)
    })
}

fn <- "fitLHS.rda"

p <- fpars <- spars <-  matrix(NA, nsim, 4)
colnames(fpars) <- colnames(spars) <- c("log.beta", "log.gamma", "log.N", "logit.i")
s.nll <- s.SSQ <- f.SSQ <- f.nll <- rep(NA, nsim)
g <- SIR.logLik()

if(file.exists(fn)){
    load(fn)
}else{
    for(i in 1:nsim){
        cat(i)
        tmp.pars <- tmpf(lhs_df[i,])
        p[i,] <- tmp.pars
        ftmp <- try(fitsir(bombay2, start = tmp.pars))
        if (!is(ftmp,"try-error")) {
            fpars[i,] <- coef(ftmp)
            f.nll[i] <- g(fpars[i,], bombay2$count)
            f.SSQ[i] <- findSSQ(bombay2, fpars[i,])$SSQ
        }
        stmp <- try(fitsir.optim(bombay2, start = tmp.pars))
        if (!is(stmp,"try-error")) {
            spars[i,] <- stmp
            s.nll[i] <- g(spars[i,], bombay2$count)
            s.SSQ[i] <- findSSQ(bombay2, spars[i,])$SSQ
        }
    }
    
    save("p", "fpars", "f.nll", "f.SSQ", "spars", "s.nll", "s.SSQ", file = fn)
}

## Another simulation

fn2 <- "fitLHS2.rda"

fpars2 <- spars2 <-  matrix(NA, nsim, 4)
colnames(fpars2) <- colnames(spars2) <- c("log.beta", "log.gamma", "log.N", "logit.i")
s.nll2 <- s.SSQ2 <- f.SSQ2 <- f.nll2 <- rep(NA, nsim)

if(file.exists(fn2)){
    load(fn2)
}else{
    for(i in 1:nsim){
        cat(i)
        ftmp <- try(fitsir(bombay2, start = fpars[i,]))
        if (!is(ftmp,"try-error")) {
            fpars2[i,] <- coef(ftmp)
            f.nll2[i] <- g(fpars2[i,], bombay2$count)
            f.SSQ2[i] <- findSSQ(bombay2, fpars2[i,])$SSQ
        }
        stmp <- try(fitsir.optim(bombay2, start = spars[i,]))
        if (!is(stmp,"try-error")) {
            spars2[i,] <- stmp
            s.nll2[i] <- g(spars2[i,], bombay2$count)
            s.SSQ2[i] <- findSSQ(bombay2, spars2[i,])$SSQ
        }
    }
    
    save("fpars2", "f.nll2", "f.SSQ2", "spars2", "s.nll2", "s.SSQ2", file = fn2)
}
