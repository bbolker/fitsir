##' Normal distribution with profiled sd
##' @param x numeric value
##' @param mean mean of distribution
##' @param log (logical) return log-likelihood
##' @return likelihood or log-likelihood vector
dnorm2 <- function(x,mean,log=FALSE) {
    rmse <- sqrt(sum((x-mean)^2)/(length(x)-1))
    return(dnorm(x,mean,sd=rmse,log=log))
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
##' @return params transformed parameter vector containing \code{beta}, \code{gamma}, \code{N}, \code{i0}
##' @export
trans.pars <- function(params) {
    tpars <- with(as.list(params),
                  c(beta=exp(log.beta),
                    gamma=exp(log.gamma),
                    N=exp(log.N),
                    i0=plogis(logit.i)))
    return(tpars)
}

##' Summarize parameters
##'
##' Generate meaningful epidemiological summary statistics (R0, r, infectious period, I(0) from a set of epidemic parameters
##' 
##' @param params parameter vector (log.beta, log.gamma, logit.i)
##' @export
##  can't call this either sum.pars or summary.pars, roxygen2 gets
##  confused ...
summarize.pars <- function(params) {
    with(as.list(params),
         c(R0=exp(log.beta-log.gamma),
           r=exp(log.beta)-exp(log.gamma),
           infper=exp(-log.gamma),
           i0=plogis(logit.i),
           I0=plogis(logit.i)*exp(log.N),
           S0=(1-plogis(logit.i))*exp(log.N),
           N=exp(log.N)))
}

summarize.pars.deriv <- function(params) {
    with(as.list(params),{
        list(
            R0.deriv=c(exp(log.beta-log.gamma), -exp(log.beta-log.gamma), 0, 0),
            r.deriv=c(exp(log.beta), -exp(log.gamma), 0, 0),
            infeper.deriv=c(0, -exp(-log.gamma), 0, 0),
            i0.deriv=c(0, 0, 0, plogis(logit.i)^2*exp(-logit.i)),
            I0.deriv=c(0, 0, plogis(logit.i)*exp(log.N), plogis(logit.i)^2*exp(-logit.i)*exp(log.N)),
            S0.deriv=c(0, 0, (1-plogis(logit.i))*exp(log.N), -plogis(logit.i)^2*exp(-logit.i)*exp(log.N)),
            N.deriv=c(0, 0, exp(log.N), 0)
        )
    })
}

##' deterministic trajectory of SIR
##' @importFrom deSolve ode
##' @param t time vector
##' @param params parameter vector (beta, gamma, N, i0)
##' @param func gradient function
##' @export
##' @useDynLib fitsir initmod
##' @examples
##' pars <- c(beta=0.4,gamma=0.2,N=5000,i0=0.001)
##' times <- 0:50
##' ss.p <- SIR.detsim(times,pars)
##' ss.i <- SIR.detsim(times,pars,type="incidence")
##' ss.d <- SIR.detsim(times,pars,type="death")
##' matplot(data.frame(ss.p,ss.i,ss.d),type = "l",xlab="time",ylab="count")
##' legend(x=0,y=800,col=1:3,lty=1:3,legend=c("prevalence","incidence","death"))
##' all.equal(cumsum(c(5, ss.i[-length(ss.i)])) - cumsum(c(0, ss.d[-length(ss.d)])), ss.p)
SIR.detsim <- function(t, params,
                       type = c("prevalence", "incidence", "death"),
                       grad = FALSE){
    type <- match.arg(type)
    
    with(as.list(params),{
            
        if(type %in% c("incidence", "death")){
            t <- c(t, 2*t[length(t)] - t[length(t)-1])
        }
        
        if (grad) {
            func <- "derivs_sens"
            yini <- c(S = N*(1-i0), logI = log(N*i0),
                nu_S_b = 0, nu_S_g = 0, nu_S_N = 1-i0, nu_S_i = -N,
                nu_I_b = 0, nu_I_g = 0, nu_I_N = i0, nu_I_i = N
            )
        }else{
            func <- "derivs"
            yini <- c(S=N*(1-i0),logI=log(N*i0))
        }
        
        odesol <- as.data.frame(ode(
            y=yini,
            times=t,
            func=func,
            parms=params[1:3],
            dllname = "fitsir",
            initfunc = "initmod",
            method = "rk4",
            hini = 0.01
        ))
        
        icol <- c("logI", "nu_I_b", "nu_I_g", "nu_I_N", "nu_I_i")
        scol <- c("S", "nu_S_b", "nu_S_g", "nu_S_N", "nu_S_i")
        
        if(type == "prevalence"){
            odesol <- odesol[,which(names(odesol) %in% icol)]
        }else{
            if(type == "death"){
                ## FIXME: type = "death", grad = TRUE sometimes returns an error
                ## In log(odesol[,"S"]) : NaNs produced
                ## possibly due to numerical imprecision... something like -1e-14
                odesol[,"logI"] <- exp(odesol[,"logI"])
                odesol <- odesol[,which(names(odesol) %in% scol)] + odesol[,which(names(odesol) %in% icol)]
            }else if(type == "incidence"){
                odesol <- odesol[,which(names(odesol) %in% scol)]
            }
            odesol <- -as.data.frame(diff(as.matrix(odesol)))
            if(!grad){
                odesol <- log(unlist(odesol, use.names = FALSE))
            }else{
                odesol[,"S"] <- log(odesol[,"S"])
            }
        }
        
        if (grad) {
            if(!all(names(odesol) == icol))
                names(odesol) <- icol
            
            logSens <- c(beta, gamma, N, i0^2*exp(-qlogis(i0)))
            odesol[,-1] <- sweep(odesol[,-1], 2, logSens, "*")
            return(odesol)
        }else{
            return(exp(odesol))
        }
    })
}

##' Log-likelihood for SIR trajectory
##' 
##' @param params parameter vector (log.N0, logit.i0, log.beta, log.gamma)
##' @param count data (epidemic counts for each time period)
##' @param times time vector
##' @param dist conditional distribution of reported data
##' @param type type of reported data
##' @param debug print debugging output?
##' @export 
SIR.logLik  <- function(params, count, times=NULL,
                  dist=c("norm", "pois", "qpois", "nbinom"),
                  type = c("prevalence", "incidence", "death"),
                  debug=FALSE) {
    dist <- match.arg(dist)
    type <- match.arg(type)
    ## HACK: nloptr appears to strip names from parameters
    ## if (is.null(params)) return(NA_real_) ## why ???
    ## if (is.null(names(params)) &&
    ##     !is.null(params)) ## ??? why ???
    ##     names(params) <- parnames(SIR.logLik)
    if (is.null(times)) times <- seq(length(count))
    if (debug) cat(params)
    tpars <- trans.pars(params)
    i.hat <- SIR.detsim(times,tpars,type)
    
    if (dist == "nbinom") {
        size <- mle.size(count,i.hat)
    } else {
        size <- NULL
    }
    r <- minusloglfun(count,i.hat,size,dist)
    if (debug) cat(" ",r,"\n")
    return(r)
}
    
## parnames() specification required in order to use
## functions with parameters specified as a
##  vector (rather than a list) with mle2

##' fitting function
##' @param data data frame
##' @param start starting parameters
##' @param dist conditional distribution of reported data
##' @param type type of reported data
##' @param method optimization method
##' @param control  control parameters for optimization
##' @param tcol column name for time variable
##' @param icol column name for count variable
##' @param debug print debugging output?
##' @export
##' @importFrom bbmle mle2
##' @examples
##' bombay2 <- setNames(bombay,c("times","count"))
##' (f1 <- fitsir(bombay2, type="death"))
##' plot(f1)
##' ss <- SIR.detsim(bombay2$times,trans.pars(coef(f1)))
##' cc <- bombay2$count
##' 
##' ## CRUDE R^2 analogue (don't trust it too far! only works if obs always>0)
##' mean((1-ss/cc)^2)
fitsir <- function(data, start=startfun(),
                   dist=c("norm", "pois", "qpois", "nbinom"),
                   type = c("prevalence", "incidence", "death"),
                   method=c("Nelder-Mead", "BFGS", "SANN"),
                   control=list(maxit=1e5),
                   tcol = "times", icol = "count",
                   timescale=NULL, ## TODO: what is this?
                   debug=FALSE) {
    dist <- match.arg(dist)
    type <- match.arg(type)
    method <- match.arg(method)
    data <- data.frame(times = data[[tcol]], count = data[[icol]])
    dataarg <- c(data,list(debug=debug, type = type, dist = dist))
    
    if (method=="BFGS" & dist=="nbinom") {
        message("Sensitivity equations are not available for negative binomial distribution.
                \nReverting back to Nelder-Mead method.")
        method <- "Nelder-Mead"
    }
    
    grad <- ifelse(method %in% c("Nelder-Mead", "SANN"), FALSE, TRUE)
    
    if (grad) {
        f.env <- new.env()
        ## set initial values
        assign("oldnll",NULL,f.env)
        assign("oldpar",NULL,f.env)
        assign("oldgrad",NULL,f.env)
        assign("data", data, f.env)
        objfun <- function(par, count, times, dist, type, debug) {
            if (identical(par,oldpar)) {
                if (debug) cat("returning old version of value\n")
                return(oldnll)
            }
            if (debug) cat("computing new version (nll)\n")
            
            v <- SIR.sensitivity(par, count, times, dist, type, debug)
            oldnll <<- v[1]
            oldgrad <<- v[-1]
            oldpar <<- par
            
            return(oldnll)
        }
        environment(objfun) <- f.env
        gradfun <- function(par, count, times, dist, type, debug) {
            if (identical(par,oldpar)) {
                if (debug) cat("returning old version of grad\n")
                return(oldgrad)
            }
            if (debug) cat("computing new version (grad)\n")
            v <- SIR.sensitivity(par, count, times, dist, type, debug)
            oldnll <<- v[1]
            oldgrad <<- v[-1]
            oldpar <<- par
            return(oldgrad)
        }
        environment(gradfun) <- f.env
        
    }else{
        objfun <- SIR.logLik
        gradfun <- NULL
    }
    attr(objfun, "parnames") <- c("log.beta","log.gamma","log.N","logit.i")
    
    m <- mle2(objfun,
              vecpar=TRUE,
              start=start,
              method=method,
              control=control,
              gr=gradfun,
              data=dataarg)
    
    m <- new("fitsir", m)
    
    if (dist == "qpois") {
        mean <- SIR.detsim(data$times, trans.pars(coef(m)), type = type)
        x <- data$count
        ss <- sum((x-mean)^2/mean)
        m@vcov <- ss/(length(x)-1) * m@vcov
    }
    
    return(m)
}

## Introducing sensitivity equations

##' Gradient of negative log likelihood with respect to each parameters
##' 
##' @param params parameter vector (log.N0, logit.i0, log.beta, log.gamma)
##' @param count data (epidemic counts for each time period)
##' @param times time vector
##' @param dist conditional distribution of reported data
##' @param type type of reported data
SIR.sensitivity <- function(params, count, times=NULL,
                     dist = c("norm", "pois", "qpois", "nbinom"),
                     type = c("prevalence", "incidence", "death"), 
                     debug = FALSE) {
    dist <- match.arg(dist)
    type <- match.arg(type)
    if (is.null(times)) times <- seq(length(count))
    tpars <- trans.pars(params)
    r <- SIR.detsim(times, tpars, type, grad = TRUE)
    
    with(as.list(c(tpars, r)), {
        i.hat <- exp(logI)
        
        nu.list <- list(nu_I_b, nu_I_g, nu_I_N, nu_I_i)
        
        nll <- minusloglfun(count, i.hat, size, dist)
        
        sensitivity <- c(nll, sapply(nu.list, function(nu) derivfun(count, i.hat, size = size, nu, dist = dist)))
        
        return(sensitivity)
    })
}

##' Derivative of negative log likelihood with respect to parameters
derivfun <- function(x,mean,size=NULL,nu,dist){
    if (dist == "qpois") dist <- "pois"
    switch(dist,
        norm = {
            n <- length(x)
            sigma2 <- sum((mean-x)^2)/(n-1)
            sum((mean-x)/sigma2 * nu + (1/(2*sigma2) - ((mean - x)^2)/(2*sigma2^2)) * sum(2 * (mean - x)/(n-1) * nu))}, 
        pois = sum((1 - x/mean) * nu)
    )
}

##' Negative log likelihood
minusloglfun <- function(x,mean,size=NULL,dist){
    if (dist == "qpois") dist <- "pois"
    switch(dist,
        norm = -sum(dnorm2(x,mean,log=TRUE)),
        pois = -sum(dpois(x,mean,log=TRUE)),
        nbinom = -sum(dnbinom(x,mu=mean,size=size,log=TRUE))
    )
}

##' Maximum likelihood estimate of negative binomial dispersion parameter
mle.size <- function(x, mean){
    nb <- function(x, mean, size) minusloglfun(x, mean, size, dist="nbinom")
    sol <- try(optim(
        par=list(size=1e2),
        fn=nb,
        gr=grad.size,
        x=x, mean=mean,
        method="BFGS",
        control=list(maxit=1e5)
    ))
    
    sol$par
}

##' derivative of nbinom nll with respect to its dispersion parameter
grad.size <- function(x, mean, size){
    n <- length(x)
    p <- size/(size + mean)
    -sum(digamma(x+size) - digamma(size) - x/(size+mean) + log(size) + 1 - log(size+mean) - size/(size+mean))
}
