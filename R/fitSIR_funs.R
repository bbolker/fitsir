
##' Starting function
##' @rdname trans.pars
##' @param log.beta log of per capita transmission rate
##' @param log.gamma log of recovery/removal rate
##' @param log.N log of population size
##' @param logit.i logit of initial proportion infectious
##' @export
startfun <- function(log.beta=log(0.12),log.gamma=log(0.09),
         log.N=log(10000),logit.i=qlogis(0.01)) {
    list(log.beta=log.beta,log.gamma=log.gamma,log.N=log.N,
         logit.i=logit.i)
}

##' Normal distribution with profiled sd
##' @param x numeric value
##' @param mean mean of distribution
##' @param log (logical) return log-likelihood?
##' @return likelihood or log-likelihood vector
dnorm2 <- function(x,mean,log=FALSE) {
    rmse <- sqrt(sum((x-mean))^2/(length(x)-1))
    return(dnorm(x,mean,sd=rmse,log=log))
}

##' gradient function (solves with {S,log(I)} for stability)
##' @param t time vector
##' @param x state vector (with names \code{S}, \code{logI})
##' @param params parameter vector (with names \code{beta} (scaled transmission rate), \code{N} (population size), \code{gamma} (recovery/removal rate)
##' @return gradient vector for a simple SIR model
SIR.grad <- function(t, x, params) {
    g <- with(as.list(c(x,params)),
          {
              c(-beta*exp(logI)*S/N,beta*S/N-gamma)
          })
    list(g)
}

##' additive log-ratio transformation and inverse
##' @param x value to transform (or inverse-transform)
alrinv <- function(x) {
    y <- exp(x)/(1+sum(exp(x)))
    c(y,1-sum(y))
}

##' @rdname alrinv
alr <- function(x) {
    y <- log(x/(1-sum(x[-length(x)])))
    y[-length(x)]
}

##' transform parameters log->exp or alr->raw
##' *assume* R=0/S=N at start of epidemic
##' @param params parameter vector containing \code{log.beta}, \code{log.gamma}, \code{log.N}, \code{logit.i}
##' @return params transformed parameter vector containing \code{beta}, \code{gamma}, \code{N}, \code{s0}, \code{i0}
trans.pars <- function(params) {
    tpars <- with(as.list(params),
                  c(beta=exp(log.beta),
                    gamma=exp(log.gamma),
                    N=exp(log.N),
                    ## use alrinv() instead?
                    s0=exp(log.N),
                    i0=plogis(logit.i)))
    if (is.list(params)) tpars <- as.list(tpars)
    return(tpars)
}

summary.pars <- function(params) {
    with(as.list(params),
         c(R0=exp(log.beta-log.gamma),
           r=exp(log.beta)-exp(log.gamma),
           infper=exp(-log.gamma),
           i0=plogis(logit.i)))
}

##' deterministic trajectory of SIR
##' @importFrom deSolve ode
##' @param t time vector
##' @param params parameter vector (beta, gamma, N, i0)
##' @param func gradient function
##' @export
SIR.detsim <- function(t, params, func=SIR.grad) {
    odesol <- with(as.list(params),
         ode(y=c(S=N, logI=log(N*i0)),
             times=t, func=func, parms=params))
    return(exp(odesol[,"logI"]))
    ## FIXME: if we want to return incidence instead
    ##        of prevalence, what is the match between
    ##        time periods and incidence? (Do we start at t=-1?)
}

##' Normal log-likelihood for SIR trajectory
##' 
##' @param params parameter vector (log.N0, logit.i0, log.beta, log.gamma)
##' @param count data (epidemic counts for each time period)
##' @param tvec time vector
##' @param debug print debugging output?
SIR.logLik <- function(params, count, tvec=NULL, debug=FALSE) {
    ## HACK: nloptr appears to strip names from parameters
    ## if (is.null(params)) return(NA_real_) ## why ???
    ## if (is.null(names(params)) &&
    ##     !is.null(params)) ## ??? why ???
    ##     names(params) <- parnames(SIR.logLik)
    if (is.null(tvec)) tvec <- seq(length(count))
    if (debug) cat(params)
    tpars <- trans.pars(params)
    i.hat <- SIR.detsim(tvec,tpars)
    r <- -sum(dnorm2(count,i.hat,log=TRUE))
    if (debug) cat(" ",r,"\n")
    return(r)
}
attr(SIR.logLik, "parnames") <- c("log.beta","log.gamma","log.N","logit.i")
## parnames() specification required in order to use
## functions with parameters specified as a
##  vector (rather than a list) with mle2

##' fitting function
##' @param data data frame with columns \code{tvec} and \code{count}
##' @param method optimization method
##' @param control  control parameters for optimization
##' @param start starting parameters
##' @param debug print debugging output?
##' @export
##' @importFrom bbmle mle2
fitsir <- function(data,method="Nelder-Mead",
                   control=list(maxit=1e5),
                   start=startfun(),debug=FALSE) {
    mle2(SIR.logLik,
         vecpar=TRUE,
         start=start,
         method=method,
         control=control,
         data=c(data,list(debug=debug)))
}

