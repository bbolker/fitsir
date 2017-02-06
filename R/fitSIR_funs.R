spline.fit <- function(tvec, count, ...){
  args <- list(...)
  
  with(args, {
    single_peak <- FALSE
    it <- 1
    spar <- 0.5
    while (!single_peak && it<itmax) {
      ss <- smooth.spline(tvec,log(count),spar=spar)
      dd <- predict(ss,deriv=1)$y
      ## change in sign of first derivative
      dds <- diff(sign(dd))
      spar <- if (spar<0.8) spar+0.05 else (1+spar)/2
      it <- it+1
      ncrit <- sum(dds<0)
      peakvals <- count[dds<0]
      relpeakvals <- peakvals[-1]/peakvals[1]
      single_peak <- ncrit==1 ||
        all(relpeakvals<relpeakcrit)
    }
    if (it==itmax) {
      ## try harder?
      stop("couldn't smooth enough")
    }
    
    return(ss)
  })
}

##' Starting function
##' @rdname trans.pars
##' @param log.beta log of per capita transmission rate
##' @param log.gamma log of recovery/removal rate
##' @param log.N log of population size
##' @param logit.i logit of initial proportion infectious
##' @export
startfun <- function(data = NULL,
                     log.beta=log(0.12),log.gamma=log(0.09),
                     log.N=log(10000),logit.i=qlogis(0.01),
                     incidence = FALSE,
                     itmax=100,relpeakcrit=0.1) {
    if (!is.null(data)) {
        tvec <- data$tvec
        count <- data$count
        ## for smooth.spline(log(count)) ...
        if (any(count<=0)) {
            count <- pmax(count,min(count[count>0])/2)
        }
        ## smooth data; start with smoothing par 0.5, try
        ## to increase it until there is a single critical point ...
        ## (check that second deriv is negative???)
        ss <- spline.fit(tvec, count, itmax = itmax, relpeakcrit = relpeakcrit)
        
        ss.data <- data.frame(tvec = tvec, count = exp(predict(ss)$y))
        
        ## find max value:
        ##  finding max based on data is a bit more robust
        ##   hard to distinguish between max and min with predict()
        ##    ufun <- function(x) predict(ss,x,deriv=1)$y
        ##    ss.tmax <- uniroot(ufun,range(tvec))$root
        ##
        ## ??? what do we do if smoothing until we have only a single
        ##   critical point 
        ss.tmax <- ss.data$tvec[which.max(ss.data$count)]
        ## find a point halfway between initial and max
        ##  scaling could be adjustable?
        ss.thalf <- min(tvec)+0.5*(ss.tmax-min(tvec))
        m1 <- lm(log(count)~tvec,data=subset(ss.data,tvec<=ss.thalf))
        r <- as.numeric(coef(m1)[2]) ##beta - gamma
        iniI <- ss.data$count[1] ## N * i0
        
        t.diff <- diff(tvec)
        t.diff <- c(t.diff[1], t.diff)
        
        if (incidence) {
            N <- cumsum(count)[length(tvec)]
            P <- ss.data$count/t.diff
            ss <- spline.fit(tvec, P, single_peak = FALSE, itmax = itmax, relpeakcrit = relpeakcrit)
        } ## if incidence
        
        ## curvature of spline at max
        ## using quadratic fit:
        ## t.sub <- (max(tvec) - ss.tmax)/2
        ## m4 <- lm(log(count)~poly(tvec,2,raw = TRUE), data = subset(ss.data, tvec> ss.tmax - t.sub & tvec < ss.tmax + t.sub))
        ## Qp.alt <- unname(2*coef(m4)[3])
        Qp.alt <- predict(ss,ss.tmax,deriv=2)$y
        if(Qp.alt > 0){
            stop("second derivative larger than 0")
        }
        Ip <- exp(max(predict(ss,tvec)$y))
        c <- -Qp.alt/Ip
        
        ss.int <- transform(ss.data, int = count * t.diff)
        ss.int <- ss.int[tvec<ss.tmax, ]
        
        d0 <- sum(ss.int[,3]) - iniI
        while(r - c * d0 < 0){
            ss.int <- ss.int[-nrow(ss.int),]
            d0 <- sum(ss.int[,3]) - iniI
        }
        
        if (incidence) {
            gamma <- 0.5 * (sqrt(4*c*N + r^2)-r)
            beta <- gamma + r
            N <- N
            d <- iniI/(t.diff[1]* beta * N)
            i0 <- 0.5 * (1-sqrt(1-4*d))
        } else {
            gamma <- c * Ip/(r - c * d0)
            beta <- gamma + r
            N <- beta*gamma/c
            i0 <- iniI/N
        }
        
        x <- c(
            log.beta = log(beta),
            log.gamma = log(gamma),
            log.N = log(N),
            logit.i = qlogis(i0)
        )
        
        return(x)
    } ## if (auto)
    
    c(log.beta=log.beta,log.gamma=log.gamma,log.N=log.N,
         logit.i=logit.i)
}

##' Normal distribution with profiled sd
##' @param x numeric value
##' @param mean mean of distribution
##' @param log (logical) return log-likelihood?
##' @return likelihood or log-likelihood vector
dnorm2 <- function(x,mean,log=FALSE) {
    rmse <- sqrt(sum((x-mean)^2)/(length(x)-1))
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
        I = exp(logI)
        dS = -beta*exp(logI)*S/N
        dlogI = beta*S/N-gamma
        list(c(dS,dlogI), I = I)
    })
}

SIR.grad.sens <- function(t, x, params) {
    g <- with(as.list(c(x,params)),
    {
        I = exp(logI)
        dS = -beta*S*I/N
        dlogI = beta*S/N-gamma
        
        grad_SS = - beta * I/N
        grad_SI = - beta * S/N
        grad_IS = beta*I/N
        grad_II = beta*S/N-gamma
	
        dnu_S_b = grad_SS * nu_S_b + grad_SI * nu_I_b - S*I/N
	
        dnu_S_N = grad_SS * nu_S_N + grad_SI * nu_I_N + beta*S*I/N^2
	
        dnu_S_g = grad_SS * nu_S_g + grad_SI * nu_I_g
	
        dnu_S_i = grad_SS * nu_S_i + grad_SI * nu_I_i
	
        dnu_I_b = grad_IS * nu_S_b + grad_II * nu_I_b + S*I/N
	
        dnu_I_N = grad_IS * nu_S_N + grad_II * nu_I_N - beta*S*I/N^2
	
        dnu_I_g = grad_IS * nu_S_g +  grad_II * nu_I_g - I
	
        dnu_I_i = grad_IS * nu_S_i + grad_II * nu_I_i
	
        list(c(dS, dlogI, dnu_S_b, dnu_S_g, dnu_S_N, dnu_S_i, dnu_I_b, dnu_I_g, dnu_I_N, dnu_I_i), I = I)
    })
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
           S0 = (1-plogis(logit.i))*exp(log.N),
           N=exp(log.N)))
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
SIR.detsim <- function(t, params, grad = FALSE,
                       incidence = FALSE){
    with(as.list(params),{
        if(incidence){
            t <- c(t[1] - diff(t[1:2]), t)
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
        
        returnName <- c("logI", "nu_I_b", "nu_I_g", "nu_I_N", "nu_I_i")
        
        if(incidence){
            odesol <- -as.data.frame(diff(as.matrix(odesol)))
            odesol[,"S"] <- log(odesol[,"S"])
            odesol <- odesol[,which(names(odesol) %in% c("S", "nu_S_b", "nu_S_g", "nu_S_N", "nu_S_i"))]
        }else{
            odesol <- odesol[,which(names(odesol) %in% returnName)]
        }
        
        if (grad) {
            if(!all(names(odesol) == returnName))
                names(odesol) <- returnName
            
            logSens <- c(beta, gamma, N, i0^2*exp(-qlogis(i0)))
            odesol[,-1] <- sweep(odesol[,-1], 2, logSens, "*")
            return(odesol)
        }else{
            return(exp(odesol))
        }
    })
}

##' Normal log-likelihood for SIR trajectory
##' 
##' @param params parameter vector (log.N0, logit.i0, log.beta, log.gamma)
##' @param count data (epidemic counts for each time period)
##' @param tvec time vector
##' @param dist conditional distribution of reported data (IGNORED)
##' @param debug print debugging output?
##' @export 
SIR.logLik  <- function(params, count, tvec=NULL,
                  dist=c("norm", "pois"),
                  debug=FALSE,
                  incidence=FALSE) {
    dist <- match.arg(dist)
    ## HACK: nloptr appears to strip names from parameters
    ## if (is.null(params)) return(NA_real_) ## why ???
    ## if (is.null(names(params)) &&
    ##     !is.null(params)) ## ??? why ???
    ##     names(params) <- parnames(SIR.logLik)
    if (is.null(tvec)) tvec <- seq(length(count))
    if (debug) cat(params)
    tpars <- trans.pars(params)
    i.hat <- SIR.detsim(tvec,tpars, incidence = incidence)
        
    r <- findNLL(count, i.hat, dist)
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
fitsir <- function(data, start=startfun(),
                   dist=c("norm", "pois"),
                   incidence=FALSE,
                   grad=FALSE,
                   method=NULL,
                   control=list(maxit=1e5),
                   timescale=NULL,
                   verbose = FALSE,
                   debug=FALSE) {
    dist <- match.arg(dist)
    dataarg <- c(data,list(debug=debug, incidence = incidence, dist = dist))
    
    ## TODO: set default method to Nelder-Mead
    ## and remove grad argument
    ## gradient evaluation should be selected based on the methods...
    
    if(grad){
        if(is.null(method)) method <- "BFGS"
        
        f.env <- new.env()
        ## set initial values
        assign("oldnll",NULL,f.env)
        assign("oldpar",NULL,f.env)
        assign("oldgrad",NULL,f.env)
        assign("data", data, f.env)
        objfun <- function(par, count, tvec, dist, incidence, debug) {
            if (identical(par,oldpar)) {
                if (verbose) cat("returning old version of value\n")
                return(oldnll)
            }
            if (verbose) cat("computing new version (nll)\n")
            
            v <- findSens(par, count, tvec, dist, incidence, debug)
            oldnll <<- v["nll"]
            oldgrad <<- v[-1]
            oldpar <<- par
            
            return(oldnll)
        }
        attr(objfun, "parnames") <- c("log.beta","log.gamma","log.N","logit.i")
        environment(objfun) <- f.env
        gradfun <- function(par, count, tvec, dist, incidence, debug) {
            if (identical(par,oldpar)) {
                if (verbose) cat("returning old version of grad\n")
                return(oldgrad)
            }
            if (verbose) cat("computing new version (grad)\n")
            v <- findSens(par, count, tvec, dist, incidence, debug)
            oldnll <<- v["nll"]
            oldgrad <<- v[-1]
            oldpar <<- par
            return(oldgrad)
        }
        environment(gradfun) <- f.env
        
        m <- mle2(objfun,
                  vecpar=TRUE,
                  start=start,
                  method=method,
                  control=control,
                  gr=gradfun,
                  data=dataarg)
        
    }else{
        if(is.null(method)) method <- "Nelder-Mead"
        
        m <- mle2(SIR.logLik,
                  vecpar=TRUE,
                  start=start,
                  method=method,
                  control=control,
                  data=dataarg)
    }
    
    ## FIXME: call mle2 only once
    
    m <- new("fitsir.mle2", m)
    
    return(m)
}

## Introducing sensitivity equations

##' returns
##'
## convert to NLL:
## NLL = C + n/2*log(SSQ/n)
## d(NLL)/dQ = d(NLL)/d(SSQ)*d(SSQ)/dQ = n/2*(n/SSQ)*1/n * d(SSQ)/dQ =
##     n/(2*SSQ) * d(SSQ)/dQ
## what is C?
findSens <- function(params, count, tvec=NULL,
                     dist = c("norm", "pois"),
                     incidence = FALSE, 
                     debug = FALSE) {
    dist <- match.arg(dist)
    if (is.null(tvec)) tvec <- seq(length(count))
    tpars <- trans.pars(params)
    r <- SIR.detsim(tvec, tpars, grad = TRUE, incidence = incidence)
    
    with(as.list(c(tpars, r)), {
        i.hat <- exp(logI)
        
        deriv_list <- list(nu_I_b, nu_I_g, nu_I_N, nu_I_i)
        
        nll <- findNLL(count, i.hat, dist)
        
        sensitivity <- c(nll, sapply(deriv_list, function(nu) findDeriv(count, i.hat, nu, dist = dist)))
        names(sensitivity) <- c("nll", "sens_beta", "sens_gamma", "sens_N", "sens_I0")
        
        return(sensitivity)
    })
}

findDeriv <- function(x, mean, nu, dist){
    switch(dist,
        norm = {
            n <- length(x)
            sigma2 <- sum((mean-x)^2)/(n-1)
            sum((mean-x)/sigma2 * nu + (1/(2*sigma2) - ((mean - x)^2)/(2*sigma2^2)) * sum(2 * (mean - x)/(n-1) * nu))}, 
        pois = sum((1 - x/mean) * nu)
    )
}

findNLL <- function(x, mean, dist){
    switch(dist,
        norm = -sum(dnorm2(x,mean,log=TRUE)),
        pois = -sum(dpois(x,mean,log=TRUE))
    )
}
