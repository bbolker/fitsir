##New parameterization:

##log.R0 <- log(R0-1)
##log.r <- log(r)
##log.S0 <- log(S0)
##log.Io <- log(I0)
startfun <- function(log.R0=log(0.12/0.09 - 1),log.r=log(0.12 - 0.09),
                     log.S0=log(10000),log.I0=log(10),auto=FALSE, incidence = FALSE,
                     data) {
  if (auto) {
    tvec <- data$tvec
    count <- data$count
    ## smooth data; start with smoothing par 0.5, try
    ## to increase it until there is a single critical point ...
    ## (check that second deriv is negative???)
    ncrit <- Inf
    it <- 1
    spar <- 0.5
    while (ncrit>1 && it<10) {
      ss <- smooth.spline(tvec,log(count),spar=spar)
      dd <- predict(ss,deriv=1)$y
      ncrit <- sum(diff(sign(dd))!=0)
      spar <- (1+spar)/2
    }
    if (it==10) stop("couldn't smooth enough")
    ss.data <- data.frame(tvec = tvec, count = exp(predict(ss)$y))
    ## find max value
    ss.tmax <- uniroot(function(x) predict(ss,x,deriv=1)$y,range(tvec))$root
    ## find a point halfway between initial and max
    ##  scaling could be adjustable?
    ss.thalf <- min(tvec)+0.5*(ss.tmax-min(tvec))
    m1 <- lm(log(count)~tvec,data=subset(ss.data,tvec<ss.thalf))
    r <- as.numeric(coef(m1)[2]) ##beta - gamma
    iniI <- ss.data$count[1] ## N * i0
    ## curvature of spline at max
    
    if(incidence){
      N = cumsum(count)[length(tvec)]
      t.diff <- diff(tvec)
      t.diff <- c(t.diff[1], t.diff)
      
      P <- ss.data$count/t.diff
      
      ncrit <- Inf
      it <- 1
      spar <- 0.5
      while (ncrit>1 && it<10) {
        ss <- smooth.spline(tvec,log(P),spar=spar)
        dd <- predict(ss,deriv=1)$y
        ncrit <- sum(diff(sign(dd))!=0)
        spar <- (1+spar)/2
      }
    }
    
    Qp.alt <- predict(ss,ss.tmax,deriv=2)$y
    Ip <- exp(max(predict(ss,tvec)$y))
    c <- -Qp.alt/Ip
    
    
    if(incidence){
      gamma = 0.5 * (sqrt(4*c*N + r^2)-r)
      beta = gamma + r
      N = N
      d = iniI/(t.diff[1]* beta * N)
      i0 = 0.5 * (1-sqrt(1-4*d))
    }else{
      gamma = -Qp.alt/r
      beta = gamma + r
      N = beta*gamma/c
      i0 = iniI/N
    }
    
    x <- list(
      log.R0 = log(beta/gamma - 1),
      log.r = log(beta-gamma),
      log.S0 = log(N*(1-i0)),
      log.I0 = log(N*i0)
    )
    
    return(x)
  }
  list(log.R0=log.R0,log.r=log.r,log.S0=log.S0,
       log.I0=log.I0)
}



dnorm2 <- function(x,mean,log=FALSE) {
  rmse <- sqrt(sum((x-mean))^2/(length(x)-1))
  return(dnorm(x,mean,sd=rmse,log=log))
}

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

trans.pars <- function(params) {
  tpars <- with(as.list(params),
                c(beta=exp(log.r - log.R0) + exp(log.r),
                  gamma=exp(log.r - log.R0),
                  N=exp(log.S0) + exp(log.I0),
                  S0=exp(log.S0),
                  I0=exp(log.I0)))
  if (is.list(params)) tpars <- as.list(tpars)
  return(tpars)
}

SIR.detsim <- function(t, params, findSens = FALSE,
                       incidence = FALSE, reportAll = FALSE){
  with(as.list(params),{
    if(incidence){
      l <- length(t)
      t.d <- diff(t[(l-1):l])
      t <- c(t, t[l]+t.d)
      t.interval <- diff(t)
    }
    
    if(findSens){
      func <- SIR.grad.sens
      yini <- c(S = S0, logI = log(I0),
                nu_beta_S = 0, nu_gamma_S = 0, nu_N_S = 1-i0,
                nu_I0_S = -N,
                nu_beta_I = 0, nu_gamma_I = 0, nu_N_I = i0,
                nu_I0_I = N)
    }else{
      func <- SIR.grad
      yini <- c(S=S0,logI=log(I0))
    }
    
    odesol <- as.data.frame(ode(y=yini,
                                times=t,
                                func=func,
                                parms=params,
                                dllname = "fitsir",
                                initfunc = "initmod"))
    
    if(findSens){
      sensName = c("nu_beta_S", "nu_gamma_S", "nu_N_S", "nu_I0_S", "nu_beta_I", "nu_gamma_I", "nu_N_I", "nu_I0_I")
      
      logSens <- c(beta, gamma, N, i0^2*exp(-qlogis(i0)))
      
      odesol[,sensName] = sweep(odesol[,sensName], 2, rep(logSens,2), "*")
    }
    
    if(reportAll){
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
      return(odesol[,"I"])
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
  attr(g, "parnames") <- c("log.R0","log.r","log.S0","log.I0")
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
  mle2(g,
       vecpar=TRUE,
       start=start,
       method=method,
       control=control,
       data=c(data,list(debug=debug)))
}

## Introducing sensitivity equations

##' Calculate SSQ
##' 
##' @param observed data and parameters for deterministic simulation
findSSQ <- function(data, params, incidence = FALSE){
  ssqL <- list()
  ssqL <- within(ssqL, {
    t = data$tvec
    sim <- SIR.detsim(t, trans.pars(params),
                      findSens = TRUE, incidence = incidence)
    obs = data$count
    pred = sim$I
    SSQ = sum(c(pred - obs)^2)
  })
  return(ssqL)
}

findSens <- function(data, params, plot.it = FALSE, log = "xy", incidence = FALSE) {
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
  return(sensitivity)
}

fitsir.optim <- function(data,
                         start = startfun(),
                         incidence = FALSE,
                         verbose = FALSE,
                         plot.it = FALSE){
  
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
                 method = "L-BFGS-B",
                 gr = gradfun)$par
  
  return(fit.p)
}

## Example

library("deSolve")
library("bbmle")
library("fitsir")

bombay2 <- setNames(bombay, c("tvec", "count"))

pars.zero <- c(
  log.R0 = log(1.5-1),
  log.r = log(0.2),
  log.S0 = log(100000),
  log.I0 = log(10)
)

start.pars <- startfun(data=bombay2, auto = TRUE)

plot(bombay2, log = "y")
lines(SIR.detsim(bombay2$tvec, trans.pars(start.pars)), col = 2)

(fpars <- coef(fitsir(data = bombay2, start = pars.zero)))
f.traj <- SIR.detsim(bombay2$tvec, trans.pars(fpars))

(fpars2 <- coef(fitsir(data = bombay2, start = startfun(data= bombay2, auto = TRUE))))
f.traj2 <- SIR.detsim(bombay2$tvec, trans.pars(fpars2))

(fpars3 <- coef(fitsir(data = bombay2, start = fpars)))
f.traj3 <- SIR.detsim(bombay2$tvec, trans.pars(fpars3))


matplot(cbind(f.traj, f.traj2, f.traj3), type = "l", log = "y")
points(bombay2, col = 4, pch = 22)
