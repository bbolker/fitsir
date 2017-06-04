##' initialize the hessian model
##' @export
initialize_hessian <- function(dist=c("gaussian", "poisson", "quasipoisson", "nbinom", "nbinom1")) {
    dist <- match.arg(dist)
    if (dist=="quasipoisson") dist <- "poisson"
    
    model <- select_model(dist)
    parnames <- c(.fitsir.pars, model@par[model@par != "param"])
    
    mu_transforms <- list(mean~.sensfun2(beta, gamma, N, i0, mean))
    mu_transforms <- gsub("mean", model@mean, as.character(mu_transforms[[1]]))
    mu_transforms <- as.formula(paste0(mu_transforms[[2]], mu_transforms[[1]], mu_transforms[[3]]))
    
    model_base <- get(paste0("loglik_", model@name, "_base"))
    
    model_hessian <- Transform(
        model_base,
        transforms=list(mu_transforms),
        par=c("beta", "gamma", "N", "i0")
    )
    
    model_hessian <- Transform(
        model_hessian,
        transforms=list(beta~exp(log.beta),
                        gamma~exp(log.gamma),
                        N~exp(log.N),
                        i0~(1+tanh(logit.i/2))/2),
        par=parnames,
        keep_hessian=TRUE
    )
    
    model_hessian
}

##' integrate second order sensitivities
##' @param t time vector
##' @param params parameter vector
SIR.detsim.hessian <- function(t, params){
    with(as.list(c(params)),{
        yini <- c(S = N*(1-i0), logI = log(N*i0),
                  nu_S_b = 0, nu_S_bb = 0, nu_S_bg = 0, nu_S_bN = 0, nu_S_bi = 0,
                  nu_S_g = 0, nu_S_gg = 0, nu_S_gN = 0, nu_S_gi = 0,
                  nu_S_N = 1-i0, nu_S_NN = 0, nu_S_Ni = -1,
                  nu_S_i = -N, nu_S_ii = 0,
                  nu_I_b = 0, nu_I_bb = 0, nu_I_bg = 0, nu_I_bN = 0, nu_I_bi = 0,
                  nu_I_g = 0, nu_I_gg = 0, nu_I_gN = 0, nu_I_gi = 0,
                  nu_I_N = i0, nu_I_NN = 0, nu_I_Ni = 1,
                  nu_I_i = N, nu_I_ii = 0
                  )
        odesol <- as.data.frame(ode(y=yini,
                                    times=t,
                                    func=SIR.grad.hessian,
                                    parms=params,
                                    method = "rk4",
                                    hini = 0.01))
    })
    
    ## TODO: allow it to work for incidence and prevalence
}

##' find Hessian
##' @param data data frame with tvec/count
##' @param params parameter vector
##' @export
SIR.hessian <- function(data, params, 
                     dist=c("gaussian", "poisson", "quasipoisson", "nbinom", "nbinom1"),
                     tcol = "times", icol = "count") {
    dist <- match.arg(dist)
    model <- initialize_hessian(dist)
    times <- data[[tcol]]
    count <- data[[icol]]
    tpars <- trans.pars(params)
    r <- SIR.detsim.hessian(times, tpars)
    
    attach(as.list(r))
    attach(as.list(params))
    
    hess <- hessian(model, count=count, mean=exp(r$logI), par=as.list(params))
    
    n <- length(params)
    
    m <- matrix(0, n, n)
    
    for(i in 1:n) {
        m[i,] <- -colSums(hess[,,i])
    }
     
    detach(as.list(r))
    detach(as.list(params))
       
    return(m)
    
}
