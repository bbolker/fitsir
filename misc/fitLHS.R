library("fitsir")
library("devtools")
source("../R/fitSIR_funs.R")
bombay2 <- setNames(bombay, c("tvec", "count"))

nsim <- 100

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
for (i in 2:ncol(lhs_df)) {
    lhs_df[,i] <- sample(lhs_df[,i])
}

p <- fpars <- spars <-  matrix(NA, nsim, 4)

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

for(i in 1:nsim){
    cat(i)
    tmp.pars <- tmpf(lhs_df[i,])
    p[i,] <- tmp.pars
    ftmp <- try(fitsir(bombay2, start = tmp.pars))
    if (!is(ftmp,"try-error")) {
        fpars[i,] <- coef(ftmp)
    }
    stmp <- try(fitsir.optim(bombay2, start = tmp.pars))
    if (!is(stmp,"try-error")) {
        spars[i,] <- stmp
    }
}

save("p", "fpars", "spars", file = "fitLHS.rda")

load("fitLHS.rda")
