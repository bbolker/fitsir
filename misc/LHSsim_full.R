library(fitsir)
library(emdbook)
source("stochsim_funs.R")
SIR.logLik <- fitsir:::SIR.logLik
fitsir.optim <- fitsir:::fitsir.optim
g <- SIR.logLik()

nsim <- 100

epi_range <- data.frame(min = c(2, 0.5, 1e3, 10),
                        max = c(3, 1, 1e5, 50),
                        row.names = c("R0", "gamma", "N", "I0"))

ltab <- as.data.frame(
    apply(epi_range, 1,
        function(x) exp(seq(log(x[1]), log(x[2]), length = nsim))
    )
)

set.seed(101)
ltab[] <- lapply(ltab, sample)

lhsf <- function(start = startfun(),
                 range = 0.5,
                 seed = NULL,
                 n = 50){
    if(!is.null(seed)) set.seed(seed)
    
    lhs.pars <- sapply(start, function(x){
        xvec <- x * seq(1 - range, 1 + range, length.out = n)
        sample(xvec)
    })
    return(lhs.pars)
}

tmppars <- startfun()
dd <- NULL

surffun <- function(x, y, 
                   data = dd,
                   llpars = tmppars){
    tmp <- llpars
    tmp["log.beta"] <- x
    tmp["logit.i"] <- y
    g(tmp, data$count, data$tvec)
}

tmpfun <- function(i, fitfun = fitsir,
                   fitrange = 0.5,
                   surfacerange = NULL,
                   ...) {
    p <- ltab[i,]
    p2 <- c(beta=unname(p["R0"]*p["gamma"]),
            p[2:3],
            i0=unname(p["I0"]/p["N"]))
    p2 <- unlist(p2)
    true.pars <- with(as.list(p2),
                      c(log.beta=log(beta),log.gamma=log(gamma),
                        log.N=log(N),logit.i=qlogis(i0)))
    d <- simfun(pars=p2,seed=101)
    
    truestart <- lhsf(true.pars, range = fitrange, seed = 101)
    
    cat("fit...\n")
    
    truefit <- apply(truestart, 1, function(x){
        f <- fitfun(d, start = x, ...)
        c(coef(f), ll = logLik(f))
    })
    
    truedf <- as.data.frame(t(truefit))
    
    cat("fit2...\n")
    
    secondfit <- apply(truedf, 1, function(x){
        f <- fitfun(d, start = x[-5], ...)
        c(coef(f), ll = logLik(f))
    })
    
    seconddf <- as.data.frame(t(secondfit))
    
    bestfit <- seconddf[which.max(seconddf$ll),]
    bestfit <- unlist(bestfit)
    
    cat("hessian...\n")
    
    hess <- findHess(d, bestfit)
    
    if(is.null(surfacerange)){
        surf <- c(0.5, 2)
    }else{
        surf <- c(1 - surfacerange, 1 + surfacerange)
    }
    
    cat("surface...\n")
    tmppars <<- bestfit[-5]; dd <<- d
    cc <- curve3d(surffun(x, y),
        xlim = sort(bestfit["log.beta"]*surf),
        ylim = sort(bestfit["logit.i"]*surf),
        n = c(51, 51)
    )
    
    return(list(
        truepars = true.pars,
        fitted = truedf,
        fitted2 = seconddf,
        bestfit = bestfit,
        hessian = hess,
        surface = cc,
        data = d
    ))
}

fn <- "LHSsim_full.rda"

resList <- resList.grad <- vector("list", nsim)

for(i in 1:nsim){
    cat(i, "\n")
    f <- try(tmpfun(i))
    cat("f-done\n")
    f.grad <- try(tmpfun(i, fitfun = fitsir.optim, nll = TRUE))
    cat("f.grad-done\n")
    if (!is(f,"try-error")) {
        resList[[i]] <- f
    }
    if (!is(f.grad,"try-error")) {
        resList.grad[[i]] <- f.grad
    }
    
    save("resList", "resList.grad", file = fn)
}

save("resList", "resList.grad", file = fn)
