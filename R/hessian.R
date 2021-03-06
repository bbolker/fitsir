##' initialize the hessian model
##' @param dist conditional distribution of reported data
##' @export
initialize_hessian <- function(dist=c("gaussian", "poisson", "quasipoisson", "nbinom", "nbinom1")) {
    dist <- match.arg(dist)
    if (dist=="quasipoisson") dist <- "poisson"
    
    model <- select_model(dist)
    parnames <- c(.upars, model@par)
    
    mu_transforms <- list(mean~.sensfun2(beta, gamma, N, i0, mean))
    mu_transforms <- gsub("mean", model@mean, as.character(mu_transforms[[1]]))
    mu_transforms <- as.formula(paste0(mu_transforms[[2]], mu_transforms[[1]], mu_transforms[[3]]))
    
    model_hessian <- Transform(
        model,
        transforms=list(mu_transforms),
        par=c("beta", "gamma", "N", "i0")
    )
    
    model_hessian <- Transform(
        model_hessian,
        transforms=list(beta~exp(log.beta),
                        gamma~exp(log.gamma),
                        N~exp(log.N),
                        i0~1/(1+exp(-logit.i0))),
        par=parnames,
        keep_hessian=TRUE
    )
    
    model_hessian
}

##' integrate second order sensitivities
##' @param t time vector
##' @param params parameter vector (beta, gamma, N, i0)
##' @param type type of count data
SIR.detsim.hessian <- function(t, params,
                               type=c("prevalence", "incidence", "death")) {
    type <- match.arg(type)
    
    with(as.list(c(params)),{
        if(type %in% c("incidence", "death")){
            t <- c(2*t[1]-t[2], t)
        }
        
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
        nn <- names(odesol)
        icol <- nn[grepl("I", nn)]
        scol <- nn[grepl("S", nn)]
        if(type == "prevalence"){
            odesol <- odesol[,which(names(odesol) %in% icol)]
        } else {
            if(type == "death"){
                odesol[,"logI"] <- exp(odesol[,"logI"])
                odesol <- odesol[,which(names(odesol) %in% scol)] + odesol[,which(names(odesol) %in% icol)]
            }else if(type == "incidence"){
                odesol <- odesol[,which(names(odesol) %in% scol)]
            }
            odesol <- -as.data.frame(diff(as.matrix(odesol)))
            odesol[,"S"] <- log(odesol[,"S"])
        }
        
        if(!all(names(odesol) == icol))
            names(odesol) <- icol
        
        return(odesol)
    })
    
    ## TODO: allow it to work for incidence and prevalence
    
}

##' find Hessian
##' @param data data frame with tvec/count
##' @param params parameter vector
##' @param dist conditional distribution of reported data
##' @param type type of reported data
##' @param tcol column name for time variable
##' @param icol column name for count variable
##' @export
SIR.hessian <- function(data, params, 
                     dist=c("gaussian", "poisson", "quasipoisson", "nbinom", "nbinom1"),
                     type=c("prevalence", "incidence", "death"),
                     tcol = "times", icol = "count") {
    type <- match.arg(type)
    dist <- match.arg(dist)
    model <- initialize_hessian(dist)
    times <- data[[tcol]]
    count <- data[[icol]]
    tpars <- trans.pars(params)
    r <- SIR.detsim.hessian(times, tpars, type=type)
    
    args <- list(object=model, count=count,mean=exp(r$logI),
                 par=params)
    args <- append(args, as.list(r)[-1])
    
    hess <- do.call(hessian, args)
    
    n <- length(params)
    
    m <- matrix(0, n, n)
    
    for(i in 1:n) {
        m[i,] <- -colSums(hess[,,i])
    }
    
    if (dist == "quasipoisson") {
        rr <- (exp(r$logI)-count)^2/exp(r$logI)
        chisq <- sum(rr)/(length(rr)-1)
        m <- m/chisq
    }
    
    return(m)
}
