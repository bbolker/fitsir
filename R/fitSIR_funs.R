
##' Starting function
##' @rdname trans.pars
##' @param log.beta log of per capita transmission rate
##' @param log.gamma log of recovery/removal rate
##' @param log.N log of population size
##' @param logit.i logit of initial proportion infectious
##' @export
startfun <- function(log.beta=log(0.12),log.gamma=log(0.09),
         log.N=log(10000),logit.i=qlogis(0.01),auto=FALSE,
         tvec,count) {
    if (auto) {
        ss <- smooth.spline(tvec,log(count),spar=0.5)
        ## find max value
        ss.tmax <- uniroot(function(x) predict(ss,x,deriv=1)$y,range(tvec))$root
        ## find a point halfway between initial and max
        ##  scaling could be adjustable?
	ss.thalf <- min(tvec)+0.5*(ss.tmax-min(tvec))
	m1 <- lm(log(count)~tvec[tvec<ss.thalf])
	r <- as.numeric(coef(m1)[2]) ##beta - gamma
	iniI <- count[1] ## N * i0
        ## curvature of spline at max
	Qp.alt <- predict(ss,ss.tmax,deriv=2)$y
	Ip <- exp(max(predict(ss,tvec)$y))
	c <- -Qp.alt/Ip
	x <- list()
	x<- within(x,{
		i0 <- 0.001 ## hack?
		N <- iniI/i0
		beta <- 0.5 *(sqrt(4*c*N + r^2)+r)
		gamma <- beta - r
	})
	return(unlist(x))

    list(log.beta=log.beta,log.gamma=log.gamma,log.N=log.N,
         logit.i=logit.i)
}

find.iniP <- function(data){
	tvec <- data$tvec
	ss <- with(bombay2,smooth.spline(tvec,log(count),spar=0.5))
	ss.tmax <- uniroot(function(x) predict(ss,x,deriv=1)$y,c(0,40))$root
	ss.thalf <- min(tvec)+(ss.tmax-min(tvec))/2
	m1 <- lm(log(count)~tvec,data=subset(bombay2,tvec<ss.thalf))
	r=as.numeric(coef(m1)[2]) ##beta - gamma
	iniI = bombay2$count[1] ## N * i0
	Qp.alt <- predict(ss,ss.tmax,deriv=2)$y
	Ip <- exp(max(predict(ss,tvec)$y))
	c = -Qp.alt/Ip
	
	x <- list()
	x<- within(x,{
		i0 = 0.001
		N = iniI/i0
		beta = 0.5 *(sqrt(4*c*N + r^2)+r)
		gamma = beta - r
	})
	return(unlist(x))
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
##' @export
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
##' @export
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
           i0=plogis(logit.i)))
}

##' deterministic trajectory of SIR
##' @importFrom deSolve ode
##' @param t time vector
##' @param params parameter vector (beta, gamma, N, i0)
##' @param func gradient function
##' @export
##' @useDynLib fitsir initmod
##' @examples
##' pars <- c(beta=0.2,gamma=0.1,N=1000,i0=0.01)
##' tvec <- seq(0,200,by=0.01)
##' ss <- SIR.detsim(tvec,pars)
##' plot(tvec,ss,type="l",xlab="time",ylab="infected")
SIR.detsim <- function(t, params, func=SIR.grad) {
    odesol <- with(as.list(params),
                   ode(y=c(S=N,logI=log(N*i0)),
                       times=t,
                       func="derivs",
                       parms=params[1:3],
                       dllname = "fitsir",
                       initfunc = "initmod",
                       nout = 1, outnames = character(0)))
    ## FIXME: bogus extra column???
    ## odesol <- with(as.list(params),
    ## ode(y=c(S=N, logI=log(N*i0)),
    ## times=t, func=func, parms=params))
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
##' @param dist conditional distribution of reported data (IGNORED)
##' @param debug print debugging output?
SIR.logLik <- function(params, count, tvec=NULL,
                       dist=dnorm2,
                       debug=FALSE) {
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
##' @examples
##' library("bbmle") ## needed at present for coef()
##' bombay2 <- setNames(bombay,c("tvec","count"))
##' ## use default starting values
##' (f1 <- fitsir(bombay2))  ## NOT a good fit
##' ss <- SIR.detsim(bombay2$tvec,trans.pars(coef(f1)))
##' cc <- bombay2$count
##' goodcoef <- c(log.beta=2.506739,log.gamma=2.475908,
##'               log.N=14.436240,logit.i=-12.782353)
##' ss2 <- SIR.detsim(bombay2$tvec,trans.pars(goodcoef))
##' plot(count~tvec,data=bombay2)
##' lines(bombay2$tvec,ss)
##' lines(bombay2$tvec,ss2,col=2)
##' ## CRUDE R^2 analogue (don't trust it too far! only works if obs always>0)
##' mean((1-ss/cc)^2)
##' mean((1-ss2/cc)^2)
fitsir <- function(data,method="Nelder-Mead",
                   control=list(maxit=1e5),
                   timescale=NULL,
                   start=startfun(),debug=FALSE) {
    mle2(SIR.logLik,
         vecpar=TRUE,
         start=start,
         method=method,
         control=control,
         data=c(data,list(debug=debug)))
}

