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
startfun <- function(log.beta=log(0.12),log.gamma=log(0.09),
                     log.N=log(10000),logit.i=qlogis(0.01),
                     auto=FALSE, incidence = FALSE,
                     data,itmax=100,relpeakcrit=0.1) {
    if (auto) {
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
        m1 <- lm(log(count)~tvec,data=subset(ss.data,tvec<ss.thalf))
        r <- as.numeric(coef(m1)[2]) ##beta - gamma
        iniI <- ss.data$count[1] ## N * i0
        ## curvature of spline at max
        
        if (incidence) {
            N = cumsum(count)[length(tvec)]
            t.diff <- diff(tvec)
            t.diff <- c(t.diff[1], t.diff)
            
            P <- ss.data$count/t.diff
            
            ss <- spline.fit(tvec, P, single_peak = FALSE, itmax = itmax, relpeakcrit = relpeakcrit)

        } ## if incidence
        
        Qp.alt <- predict(ss,ss.tmax,deriv=2)$y
        Ip <- exp(max(predict(ss,tvec)$y))
        c <- -Qp.alt/Ip
        
        if (incidence) {
            gamma <- 0.5 * (sqrt(4*c*N + r^2)-r)
            beta <- gamma + r
            N <- N
            d <- iniI/(t.diff[1]* beta * N)
            i0 <- 0.5 * (1-sqrt(1-4*d))
        } else {
            gamma <- -Qp.alt/r
            beta <- gamma + r
            N <- beta*gamma/c
            i0 <- iniI/N
        }
        
        x <- list(
            log.beta = log(beta),
            log.gamma = log(gamma),
            log.N = log(N),
            logit.i = qlogis(i0)
        )
        
        return(x)
    } ## if (auto)
    
    list(log.beta=log.beta,log.gamma=log.gamma,log.N=log.N,
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
	
        dnu_beta_S = grad_SS * nu_beta_S + grad_SI * nu_beta_I - S*I/N
	
        dnu_N_S = grad_SS * nu_N_S + grad_SI * nu_N_I + beta*S*I/N^2
	
        dnu_gamma_S = grad_SS * nu_gamma_S + grad_SI * nu_gamma_I
	
        dnu_I0_S = grad_SS * nu_I0_S + grad_SI * nu_I0_I
	
        dnu_beta_I = grad_IS * nu_beta_S + grad_II * nu_beta_I + S*I/N
	
        dnu_N_I = grad_IS * nu_N_S + grad_II * nu_N_I - beta*S*I/N^2
	
        dnu_gamma_I = grad_IS * nu_gamma_S +  grad_II * nu_gamma_I - I
	
        dnu_I0_I = grad_IS * nu_I0_S + grad_II * nu_I0_I
	
        list(c(dS, dlogI, dnu_beta_S, dnu_gamma_S, dnu_N_S, dnu_I0_S, dnu_beta_I, dnu_gamma_I, dnu_N_I, dnu_I0_I), I = I)
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
SIR.detsim <- function(t, params, findSens = FALSE,
                       incidence = FALSE, reportAll = FALSE){
    with(as.list(params),{
        if(incidence){
            l <- length(t)
            t.d <- diff(t[(l-1):l])
            t <- c(t, t[l]+t.d)
            t.interval <- diff(t)
        }
        
        if (findSens) {
            func <- "derivs_sens"
            yini <- c(S = N*(1-i0), logI = log(N*i0),
                      nu_beta_S = 0, nu_gamma_S = 0, nu_N_S = 1-i0,
                      nu_I0_S = -N,
                      nu_beta_I = 0, nu_gamma_I = 0, nu_N_I = i0,
                      nu_I0_I = N)
        }else{
            func <- "derivs"
            yini <- c(S=N*(1-i0),logI=log(N*i0))
        }
        
        odesol <- as.data.frame(ode(y=yini,
                                    times=t,
                                    func=func,
                                    parms=params[1:3],
                                    dllname = "fitsir",
                                    initfunc = "initmod",
                                    method = "rk4",
                                    hini = 0.01))
        
        if (findSens) {
            sensName = c("nu_beta_S", "nu_gamma_S", "nu_N_S", "nu_I0_S", "nu_beta_I", "nu_gamma_I", "nu_N_I", "nu_I0_I")
            logSens <- c(beta, gamma, N, i0^2*exp(-qlogis(i0)))
            odesol[,sensName] = sweep(odesol[,sensName], 2, rep(logSens,2), "*")
        }
        
        if (reportAll){
            return(odesol)
        }
        
        if(incidence){
            odesol <- as.data.frame(diff(as.matrix(odesol)))
            if(findSens){
                odesol <- -odesol[,c("S", "nu_beta_S", "nu_gamma_S", "nu_N_S", "nu_I0_S")]
                names(odesol) <- c("I", "nu_beta_I", "nu_gamma_I", "nu_N_I", "nu_I0_I")		
            }else{
                odesol <- -odesol
                names(odesol) <- c("NA.t","I","NA.p")
            }
        }
        
        if(findSens){
            return(odesol)
        }else{
            return(exp(odesol[,"logI"]))
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
SIR.logLik <- function(incidence = FALSE){
    g <- function(params, count, tvec=NULL,
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
        i.hat <- SIR.detsim(tvec,tpars, incidence = incidence)
        
        r <- -sum(dnorm2(count,i.hat,log=TRUE))
        if (debug) cat(" ",r,"\n")
        return(r)
    }
    attr(g, "parnames") <- c("log.beta","log.gamma","log.N","logit.i")
    return(g)
}
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
                   incidence = FALSE,
                   start=startfun(),debug=FALSE) {
    g <- SIR.logLik(incidence = incidence)
    m <- mle2(g,
         vecpar=TRUE,
         start=start,
         method=method,
         control=control,
         data=c(data,list(debug=debug)))
    ## FIXME: extend S4 class?
    ## class(m) <- c("fitsir","mle2")
    return(m)
}

## Introducing sensitivity equations

##' Calculate SSQ
##' 
##' @param data data for deterministic simulation (data frame with columns \code{tvec} and \code{count})
##' @param params parameter vector
##' @param incidence (logical) compute SSQ based on incidence?
##' @param SSQonly (logical) return only the sum of squares? (alternatively, return derivatives of SSQ with respect to parameters)
findSSQ <- function(data, params, incidence = FALSE, SSQonly = FALSE){
    ssqL <- list()
    ssqL <- within(ssqL, {
        t <- data$tvec
        sim <- SIR.detsim(t, trans.pars(params),
                          findSens = TRUE, incidence = incidence)
        obs <- data$count
        pred <- exp(sim$logI)
        SSQ <- sum((pred - obs)^2)
    })
    if (SSQonly) return(ssqL$SSQ)
    return(ssqL)
}

##' returns
##'
## convert to NLL:
## NLL = C + n/2*log(SSQ/n)
## d(NLL)/dQ = d(NLL)/d(SSQ)*d(SSQ)/dQ = n/2*(n/SSQ)*1/n * d(SSQ)/dQ =
##     n/(2*SSQ) * d(SSQ)/dQ
## what is C?
findSens <- function(data, params, plot.it = FALSE, log = "xy",
                     incidence = FALSE, sensOnly = FALSE,
                     nll = FALSE) {
    ssqL <- findSSQ(data, params, incidence = incidence)
    if (plot.it) {
        matplot(cbind(ssqL$obs, ssqL$pred), log=log, type = "l")
    }
    dSSQ <- 2 * (ssqL$pred - ssqL$obs)
    sensitivity <- with(ssqL$sim,
                        c(
                            SSQ = ssqL$SSQ,
                            SSQ_beta = sum(dSSQ * nu_beta_I),
                            SSQ_gamma = sum(dSSQ * nu_gamma_I),
                            SSQ_N = sum(dSSQ * nu_N_I),
                            SSQ_I0 = sum(dSSQ * nu_I0_I)
                        ))
    if (nll) {
        n <- length(ssqL$SSQ)
        NLL <- n/2*log(sensitivity["SSQ"]/n)
        NLLsens <- sensitivity[-1]*n/(2*sensitivity["SSQ"])
        sensitivity <- c(NLL=unname(NLL),
                         setNames(NLLsens,
                                  paste("SSQ",c("beta","gamma","N","I0"),
                                         sep="_")))
    }
    if(sensOnly){
        sensitivity <- sensitivity[-1]
    }
    return(sensitivity)
}

fitsir.optim <- function(data,
                         start = startfun(),
                         incidence = FALSE,
                         verbose = FALSE,
                         plot.it = FALSE,
                         debug = FALSE){
    
    if(plot.it){
        plot(data)
    }
    
    f.env <- new.env()
    ## set initial values
    assign("oldSSQ",NULL,f.env)
    assign("oldpar",NULL,f.env)
    assign("oldgrad",NULL,f.env)
    assign("data", data, f.env)
    objfun <- function(par) {
        if (identical(par,oldpar)) {
            if (verbose) cat("returning old version of SSQ\n")
            return(oldSSQ)
        }
        if(debug) cat(par, "\n")
        if (verbose) cat("computing new version (SSQ)\n")
        
        v <- findSens(data, par, incidence = incidence)
        oldSSQ <<- v["SSQ"]
        oldgrad <<- v[-1]
        oldpar <<- par
        
        if(plot.it){
            lines(data$tvec,
                  SIR.detsim(data$tvec, trans.pars(par),
                             incidence = incidence))
        }
        return(oldSSQ)
    }
    environment(objfun) <- f.env
    gradfun <- function(par) {
        if (identical(par,oldpar)) {
            if (verbose) cat("returning old version of grad\n")
            return(oldgrad)
        }
        if (verbose) cat("computing new version (grad)\n")
        v <- findSens(data, trans.pars(par), incidence = incidence)
        oldSSQ <<- v["SSQ"]
        oldgrad <<- v[-1]
        oldpar <<- par
        return(oldgrad)
    }
    environment(gradfun) <- f.env
    
    fit.p <- optim(fn = objfun,
                   par = start,
                   method = "BFGS",
                   gr = gradfun)$par
    
    return(fit.p)
}
