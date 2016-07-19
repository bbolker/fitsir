if (file.exists("~/Rlibs")) .libPaths("~/Rlibs")

library(fitsir)
library(plyr)
source("stochsim_funs.R")

## beta_range <- c(0.16,0.75)
## sample over R0 and post-process to beta=R0/gamma, to
##   avoid R0<1 cases from sampling beta and gamma independently ...
R0_range <- c(1.1,8)
gamma_range <- c(0.07,0.3)
N_range <- c(500,3e4)
I0_range <- c(1,20)


nsim <- 1000
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
fn <- "stochsim_4.rda"
nmvec <- nmvec2 ## from stochsim_funs.R
null_out <- setNames(rep(NA,length(nmvec)),nmvec)

## could use aaply but for loop is more transparent
res <- matrix(NA,nrow=nsim,ncol=length(nmvec),
              dimnames=list(NULL,nmvec))

for (i in 1:nsim) {
    cat(i,"\n")
    p <- lhs_df[i,]
    p2 <- c(p[2:3],
            beta=unname(p["R0"]*p["gamma"]),
            i0=unname(p["I0"]/p["N"]))
    d <- simfun(pars=p2,tmax=100,dt=1,rpars=list(size=3),seed=101)
    ss0 <- startfun(auto=TRUE,data=d)
    f <- try(fitfun2(d))
    if (!is(f,"try-error")) {
        res[i,] <- f
    }
    save("res",file=fn)
}
save("res","lhs_df",file=fn)

    
