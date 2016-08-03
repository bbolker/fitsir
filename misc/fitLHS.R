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

fn <- "fitLHS.rda"

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
        }
        stmp <- try(fitsir.optim(bombay2, start = tmp.pars))
        if (!is(stmp,"try-error")) {
            spars[i,] <- stmp
        }
    }
    
    save("p", "fpars", "spars", file = fn)
}

colnames(fpars) <- colnames(spars) <- c("log.beta", "log.gamma", "log.N", "logit.i")

matplot(apply(fpars, 2, sort), type = "l")
matplot(apply(spars, 2, sort), type = "l")

##terrible fit
findSens(bombay2, fpars[5,], plot.it = TRUE)
findSens(bombay2, spars[5,], plot.it = TRUE)
##running it once more brings it to a local minimum
tmp.pars <- fitsir.optim(bombay2, start = spars[5,], plot.it = TRUE)




