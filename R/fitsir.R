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
##' Generate meaningful epidemiological summary statistics (R0, r, infectious period, I(0)) from a set of epidemic parameters
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
                ## possibly due to underflow
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
            
            logSens <- c(beta, gamma, N, i0*(1-i0))
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
                        dist=c("gaussian", "poisson", "quasipoisson", "nbinom", "nbinom1"),
                        type = c("prevalence", "incidence", "death"),
                        debug=FALSE) {
    dist <- match.arg(dist)
    model <- select_model(dist)
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
    
    if (grepl("nbinom", dist)) {
        par <- params[5]
    } else {
        par <- NULL
    }
    r <- -sum(Eval(model, count, i.hat, par))
    if (debug) cat(" ",r,"\n")
    return(r)
}

.fitsir.pars <- c("log.beta","log.gamma","log.N","logit.i")
  
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
##' @seealso startfun
##' @export
##' @importFrom bbmle mle2
##' @examples
##' bombay2 <- setNames(bombay,c("times","count"))
##' (f1 <- fitsir(bombay2, start=c(log.beta=log(2), log.gamma=log(1), log.N=log(1000), logit.i=qlogis(0.001)), type="death"))
##' plot(f1)
##' ss <- SIR.detsim(bombay2$times,trans.pars(coef(f1)))
##' cc <- bombay2$count
##' 
##' ## CRUDE R^2 analogue (don't trust it too far! only works if obs always>0)
##' mean((1-ss/cc)^2)
fitsir <- function(data, start,
                   dist=c("gaussian", "poisson", "quasipoisson", "nbinom", "nbinom1"),
                   type = c("prevalence", "incidence", "death"),
                   method=c("Nelder-Mead", "BFGS", "SANN"),
                   control=list(maxit=1e5),
                   tcol = "times", icol = "count",
                   debug=FALSE,
                   ...) {
    dist <- match.arg(dist)
    type <- match.arg(type)
    method <- match.arg(method)
    data <- data.frame(times = data[[tcol]], count = data[[icol]])
    dataarg <- c(data,list(debug=debug, type = type, dist=dist))
    
    model <- select_model(dist)
    parnames <- c(.fitsir.pars, model@par[model@par != "param"])
    
    if(missing(start) | length(start) != length(parnames)) {
        stop.msg <- paste0("'start' must be a named vector with following parameters:\n", paste(parnames, collapse = ', '))
        stop(stop.msg)
    } else if (!all(parnames %in% names(start))) {
        stop.msg2 <- paste0("names of 'start' does not match names of the parameters (", paste(parnames, collapse = ', '), ").
                           \n")
        stop(stop.msg2)
    }
    
    if (method=="BFGS") {
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
    attr(objfun, "parnames") <- parnames
    
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
    }
    
    return(m)
}

## Introducing sensitivity equations

##' Gradient of negative log likelihood with respect to each parameters
##' 
##' @param params parameter vector (log.N, logit.i, log.beta, log.gamma)
##' @param count data (epidemic counts for each time period)
##' @param times time vector
##' @param dist conditional distribution of reported data
##' @param type type of reported data
##' @param debug print debugging output?
##' @examples
##' fitsir:::SIR.sensitivity(c(log.beta=2,log.gamma=0,logit.i=-3,log.N=4),
##'                  count=c(1,2,4,7,3),
##'                  times=1:5,
##'                  dist="gaussian",
##'                  type="prevalence")
##' 
##' @export
SIR.sensitivity <- function(params, count, times=NULL,
                            dist=c("gaussian", "poisson", "quasipoisson", "nbinom", "nbinom1"),
                            type = c("prevalence", "incidence", "death"),
                            debug=FALSE) {
    dist <- match.arg(dist)
    model <- select_model(dist)
    type <- match.arg(type)
    if (is.null(times)) times <- seq(length(count))
    tpars <- trans.pars(params)
    r <- SIR.detsim(times, tpars, type, grad = TRUE)
    ## FIXME: Make this more general
    if (grepl("nbinom", dist)) {
        par <- params[5]
    } else {
        par <- NULL
    }
    
    res <- with(as.list(c(tpars, r)), {
        i.hat <- exp(logI)
        nu.list <- list(nu_I_b, nu_I_g, nu_I_N, nu_I_i)
        nll <- -sum(Eval(model, count, i.hat, par))
        sensitivity <- c(nll, sapply(nu.list, function(nu) -sum(grad(model,count,i.hat,par,param=NULL,nu=nu,var="param")[[1]])))
        names(sensitivity) <- c("value",.fitsir.pars)
        if (grepl("nbinom", dist)) {
            sensitivity <- c(sensitivity,
                             -sum(grad(model,count,i.hat,par,param=NULL,nu=NULL,var=1)[[1]]))
            names(sensitivity)[-c(1:5)] <- model@par[model@par != "param"]
        }
            
        return(sensitivity)
    })
    return(res)
}

##' Select likelihood model
##' @param dist conditional distribution of reported data
##' @export
select_model <- function(dist = c("gaussian", "poisson", "quasipoisson", "nbinom", "nbinom1")) {
    dist <- match.arg(dist)
    if (dist == "quasipoisson") dist <- "poisson"
    get(paste0("loglik_", dist))
}

##' Negative log likelihood functions
##' @param x vector of observations
##' @param mean vector of means
##' @param par additional parameter
##' @param dist conditional distribution of reported data
minusloglfun <- function(x,mean,par,dist){
    if (dist == "quasipoisson") dist <- "poisson"
    model <- get(paste0("loglik_", dist))
    -sum(Eval(model, x, mean, par))
}

##' Maximum likelihood estimate of negative binomial dispersion parameter
mledsp <- function(x,mean,dist=c("nbinom", "nbinom1")){
    dist <- match.arg(dist)
    model <- get(paste0("loglik_", dist))
    
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
