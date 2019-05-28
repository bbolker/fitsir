##' transform parameters
##' 
##' @param params parameter vector containing \code{log.beta}, \code{log.gamma}, \code{log.N}, \code{logit.i}
##' @return params transformed parameter vector containing \code{beta}, \code{gamma}, \code{N}, \code{i0}
##' @export
trans.pars <- function(params, scale=c("constrained", "unconstrained")) {
    scale <- match.arg(scale)
    tpars <- switch(scale,
        constrained=with(as.list(params),
                         c(beta=exp(log.beta),
                           gamma=exp(log.gamma),
                           N=exp(log.N),
                           i0=plogis(logit.i0))),
        unconstrained=with(as.list(params),
                           c(log.beta=log(beta),
                             log.gamma=log(gamma),
                             log.N=log(N),
                             logit.i0=qlogis(i0))))
    tpars
}

trans.pars.loglik <- function(params, model, scale=c("constrained", "unconstrained")) {
    scale <- match.arg(scale)
    tpars <- switch(scale,
        constrained=with(as.list(params), {
            switch(model@name,
                   gaussian=c(sigma=exp(log.sigma)),
                   nbinom=c(k=exp(log.k)),
                   nbinom1=c(phi=exp(log.phi))
            )
        }),
        unconstrained=with(as.list(params), {
            switch(model@name,
                   gaussian=c(log.sigma=log(sigma)),
                   nbinom=c(log.k=log(k)),
                   nbinom1=c(log.phi=log(phi))
            )
        }))
    tpars
}

.linkfun <- function(params, model) {
    tpars <- c(trans.pars(params, "constrained"), 
      trans.pars.loglik(params, model, "constrained"))
    
    mu.eta <- with(as.list(params),{
        c(log.beta=exp(log.beta),
          log.gamma=exp(log.gamma),
          log.N=exp(log.N),
          logit.i0=.logit.eta(logit.i0))
    })
    
    mu.eta.loglik <- with(as.list(params), {
        switch(model@name,
               gaussian=c(log.sigma=exp(log.sigma)),
               nbinom=c(log.k=exp(log.k)),
               nbinom1=c(log.phi=exp(log.phi))
        )
    })
    
    list(
        param=tpars,
        mu.eta=c(mu.eta, mu.eta.loglik)
    )
}

.logit.eta <- make.link("logit")$mu.eta

##' Summarize parameters
##'
##' Generate meaningful epidemiological summary statistics (R0, r, infectious period, I(0)) from a set of epidemic parameters
##' 
##' @param params parameter vector (log.beta, log.gamma, logit.i)
##' @export summarize.pars
##  can't call this either sum.pars or summary.pars, roxygen2 gets
##  confused ...
summarize.pars <- function(params) {
    with(as.list(params),
         c(R0=exp(log.beta-log.gamma),
           r=exp(log.beta)-exp(log.gamma),
           infper=exp(-log.gamma),
           i0=plogis(logit.i0),
           I0=plogis(logit.i0)*exp(log.N),
           S0=(1-plogis(logit.i0))*exp(log.N),
           N=exp(log.N)))
}

##' @rdname summarize.pars
summarize.pars.jacobian <- function(params) {
    with(as.list(params),{
        m <- list(
            R0.deriv=c(exp(log.beta-log.gamma), -exp(log.beta-log.gamma), 0, 0),
            r.deriv=c(exp(log.beta), -exp(log.gamma), 0, 0),
            infeper.deriv=c(0, -exp(-log.gamma), 0, 0),
            i0.deriv=c(0, 0, 0, plogis(logit.i)^2*exp(-logit.i)),
            I0.deriv=c(0, 0, plogis(logit.i)*exp(log.N), plogis(logit.i)^2*exp(-logit.i)*exp(log.N)),
            S0.deriv=c(0, 0, (1-plogis(logit.i))*exp(log.N), -plogis(logit.i)^2*exp(-logit.i)*exp(log.N)),
            N.deriv=c(0, 0, exp(log.N), 0)
        )
        unname(do.call(rbind, m))
    })
}

##' deterministic trajectory of SIR
##' @importFrom deSolve ode
##' @param t time vector
##' @param params parameter vector (beta, gamma, N, i0)
##' @param type type of count data
##' @param grad (logical) return gradient with respect to unconstrained parameters
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
            t <- c(2*t[1]-t[2], t)
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
            method = "rk4"
        ))
        
        icol <- c("logI", "nu_I_b", "nu_I_g", "nu_I_N", "nu_I_i")
        scol <- c("S", "nu_S_b", "nu_S_g", "nu_S_N", "nu_S_i")
        
        if(type == "prevalence"){
            odesol <- odesol[,which(names(odesol) %in% icol)]
        }else{
            if(type == "death"){
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
            names(odesol) <- c("logI", "beta", "gamma", "N", "i0")
            
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
##' @param model log likelihood model
##' @param type type of reported data
##' @param debug print debugging output?
##' @export 
SIR.logLik  <- function(params, count, times=NULL,
                        model,
                        type = c("prevalence", "incidence", "death"),
                        debug=FALSE) {
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
    
    r <- -sum(Eval(model, count, i.hat, params[-c(1:4)]))
    
    if (debug) cat(" ",r,"\n")
    return(r)
}

##' constrained parameters
.cpars <- c("beta", "gamma", "N", "i0")

##' unconstrained parameters
.upars <- c("log.beta","log.gamma","log.N","logit.i0")

  
##' fitting function
##' @param data data frame
##' @param start starting parameters
##' @param dist conditional distribution of reported data
##' @param type type of reported data
##' @param method optimization method (see \link{mle2})
##' @param optimizer optimizer to be used (see \link{mle2})
##' @param control  control parameters for optimization
##' @param tcol column name for time variable
##' @param icol column name for count variable
##' @param debug print debugging output?
##' @param ... Further arguments to pass to optimizer
##' @seealso startfun mle2
##' @export
##' @importFrom bbmle mle2
##' @examples
##' harbin2 <- setNames(harbin,c("times","count"))
##' (f1 <- fitsir(harbin2, 
##'               start=c(beta=2, gamma=1, N=2e3, i0=0.0001, sigma=10),
##'               type="death"))
##' plot(f1)
##' 
##' ## CRUDE R^2 analogue (don't trust it too far!)
##' ss <- SIR.detsim(harbin2$times,trans.pars(coef(f1)), type="death")
##' cc <- harbin2$count
##' 
##' cor(ss, cc)^2
##'
##' f1_g2 <- fitsir(harbin2, 
##'              start=c(beta=2, gamma=1, N=2e3, i0=0.0001),
##'               family="gaussian2",
##'               type="death")
##' plot(f1)
fitsir <- function(data, start,
                   dist=c("gaussian2",
                          "gaussian", "poisson", "quasipoisson",
                          "nbinom", "nbinom1"),
                   type = c("prevalence", "incidence", "death"),
                   method="BFGS",
                   control=list(maxit=1e5),
                   tcol = "times", icol = "count",
                   debug=FALSE,
                   ...) {
    dist <- match.arg(dist)
    type <- match.arg(type)
    data <- data.frame(times = data[[tcol]], count = data[[icol]])
    
    model <- select_model(dist)
    
    if(method=="mle") {
        
    }
    
    parnames <- c(.cpars, model@par)
    
    if(missing(start) | length(start) != length(parnames)) {
        stop.msg <- paste0("'start' must be a named vector with following parameters:\n", paste(parnames, collapse = ', '))
        stop(stop.msg)
    } else if (!all(parnames %in% names(start))) {
        stop.msg2 <- paste0("names of 'start' does not match names of the parameters (", paste(parnames, collapse = ', '), ").
                            \n")
        stop(stop.msg2)
    }
    
    ## put parameters in the right order
    sir.start <- start[.cpars]
    loglik.start <- start[which(!(names(start) %in% .cpars))]
    
    start <- c(
        trans.pars(sir.start, "unconstrained"),
        trans.pars.loglik(loglik.start, model, "unconstrained")
    )
    
    dataarg <- c(data,list(debug=debug, type = type, model=model))
    
    f.env <- new.env()
    ## set initial values
    assign("oldnll",NULL,f.env)
    assign("oldpar",NULL,f.env)
    assign("oldgrad",NULL,f.env)
    assign("data", data, f.env)
    objfun <- function(par, count, times, model, type, debug) {
        if (identical(par,oldpar)) {
            if (debug) cat("returning old version of value\n")
            return(oldnll)
        }
        if (debug) cat("computing new version (nll)\n")
        
        v <- SIR.sensitivity(par, count, times, model, type, debug)
        oldnll <<- v[1]
        oldgrad <<- v[-1]
        oldpar <<- par
        
        return(oldnll)
    }
    environment(objfun) <- f.env
    gradfun <- function(par, count, times, model, type, debug) {
        if (identical(par,oldpar)) {
            if (debug) cat("returning old version of grad\n")
            return(oldgrad)
        }
        if (debug) cat("computing new version (grad)\n")
        v <- SIR.sensitivity(par, count, times, model, type, debug)
        oldnll <<- v[1]
        oldgrad <<- v[-1]
        oldpar <<- par
        return(oldgrad)
    }
    environment(gradfun) <- f.env
    
    attr(objfun, "parnames") <- names(start)
    
    m <- mle2(objfun,
              vecpar=TRUE,
              start=start,
              method=method,
              control=control,
              gr=gradfun,
              data=dataarg,
              ...)
    
    m <- new("fitsir", m)
    
    if (dist == "quasipoisson") {
        rr <- residuals(m)
        chisq <- sum(rr)/(length(rr)-1)
        m@vcov <- chisq * m@vcov
        m@details$hessian <- m@details$hessian/chisq
    }
    
    return(m)
}

## Introducing sensitivity equations

##' Gradient of negative log likelihood with respect to each parameters
##' 
##' @param params parameter vector (log.N, logit.i, log.beta, log.gamma)
##' @param count data (epidemic counts for each time period)
##' @param times time vector
##' @param model log-likelihood model
##' @param type type of reported data
##' @param debug print debugging output?
SIR.sensitivity <- function(params, count, times=NULL,
                            model,
                            type = c("prevalence", "incidence", "death"),
                            debug=FALSE) {
    type <- match.arg(type)
    if (is.null(times)) times <- seq(length(count))
    linkpar <- .linkfun(params, model)
    
    tpars <- linkpar$param
    r <- SIR.detsim(times, tpars, type, grad = TRUE)
    
    mean <- exp(r$logI)
    
    loglik.par <- as.list(tpars[-c(1:4)])
    
    nll <- -sum(Eval(model, count, mean, loglik.par))
    loglik.gr <- grad(model, count, mean, loglik.par, var=c(model@mean, model@par))
    
    sensitivity <- -colSums(loglik.gr[[1]] * r[,-1])
    
    if(length(loglik.gr) > 1) sensitivity <- c(sensitivity, -sapply(loglik.gr[-1], sum))
    
    sensitivity <- sensitivity * linkpar$mu.eta
    
    c(nll, sensitivity)
}

##' Maximum likelihood estimate of negative binomial dispersion parameter
##' @param x vector of counts
##' @param mean mean of the distribution
##' @param dist conditional distribution of reported data
mledsp <- function(x,mean,dist=c("nbinom", "nbinom1")){
    dist <- match.arg(dist)
    model <- select_model(dist)
    
    start <- switch(dist,
        "nbinom"=log(1e2),
        "nbinom1"=log(2)
    )
    
    fn <- function(x, mean, par) -sum(Eval(model,x,mean, par))
    
    gr <- function(x, mean, par) -sum(grad(model,x,mean,par,param=NULL,nu=NULL,var=1)[[1]])
    
    sol <- try(optim(
        par=list(par=start),
        fn=fn,
        gr=gr,
        x=x, mean=pmax(mean, 1e-100),
        method="BFGS",
        control=list(maxit=1e4)
    ))
    
    if(inherits(sol, "try-error")) {
        browser()
    }
    sol$par
}
