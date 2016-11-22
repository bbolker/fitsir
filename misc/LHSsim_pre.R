library(fitsir)
library(emdbook)
library(rgl)
source("stochsim_funs.R")

truepars <- c(log.beta = log(2),
    log.gamma = log(0.8),
    log.N = log(1e6),
    logit.i = qlogis(1e-5)
)

dd <- simfun(trans.pars(truepars), seed = 123, rpars = list(size = 4))

lhsf <- function(start = truepars,
                 range = 0.2,
                 seed = NULL,
                 n = 50){
    if(!is.null(seed)) set.seed(seed)
    
    lhs.pars <- sapply(start, function(x){
        xvec <- x * seq(1 - range, 1 + range, length.out = n)
        sample(xvec)
    })
    return(lhs.pars)
}

lhspars <- lhsf(seed = 101)

fn <- "LHSsim_pre.rda"

## fitting... starting values are sampled using LHS

if(file.exists(fn)){
    load(fn)
}elsle{
    resList <- apply(lhspars, 1, function(x){
        fitsir(dd, start = x)
    })
    save("lhspars", "resList", file = fn)
}

fitdf <- do.call("rbind", lapply(resList, function(x){
    c(coef(x), ll = logLik(x))
}))

## it seems like most parameters converge to one place
## one of them is getting stuck?
## maybe try this with auto start or with wider range?
pairs(fitdf)

SIR.logLik <- fitsir:::SIR.logLik()

## Assuming that we know N and gamma...
tmpfun <- function(x, y, 
                   data = dd,
                   llpars = truepars){
    llpars["log.beta"] <- x
    llpars["logit.i"] <- y
    SIR.logLik(llpars, data$count, data$tvec)
}

cc <- curve3d(tmpfun(x, y),
    xlim = truepars["log.beta"]*c(0.4, 1.6),
    ylim = truepars["logit.i"]*c(1.6, 0.4),
    n = c(51, 51)
)


image(cc, xlab = "log.beta", ylab = "logit.i")
persp3d(cc, col = "blue")

## might be better to use fitted parameters
bestfit <- fitdf[which.max(fitdf[,5]),-5]

cc2 <- curve3d(tmpfun(x, y, llpars = bestfit, data = dd),
    xlim = bestfit["log.beta"]*c(0.4, 1.6),
    ylim = bestfit["logit.i"]*c(1.6, 0.4),
    n = c(51, 51)
)

image(cc2, xlab = "log.beta", ylab = "logit.i")
persp3d(cc2, col = "blue")

## trying with bombay data

bombay2 <- setNames(bombay, c("tvec", "count"))

fitb <- fitsir(bombay2)
fpars <- coef(fitb)

cc_bombay <- curve3d(tmpfun(x, y, llpars = fpars, data = bombay2),
    xlim = fpars["log.beta"]*c(0.9, 1.1),
    ylim = fpars["logit.i"]*c(1.1, 0.9),
    n = c(81, 81)
)

image(cc_bombay, xlab = "log.beta", ylab = "logit.i")
persp3d(cc_bombay, col = "blue")


save(list=ls(pattern="^cc*"),file="LHSsim_pre_cc.rda")
