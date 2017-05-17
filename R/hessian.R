initialize_hessian <- function(dist=c("gaussian", "poisson", "quasipoisson", "nbinom", "nbinom1")) {
    dist <- match.arg(dist)
    if (dist=="quasipoisson") dist <- "poisson"
    model <- switch (dist,
        gaussian={
            Transform(
                loglik_gaussian_base,
                transforms = list(mu~.sensfun2(param1, param2, mu)),
                par=c("param1", "param2")
            )
            
        },
        poisson={
            Transform(
                loglik_poisson_base,
                transforms = list(lambda~.sensfun2(param1, param2, lambda)),
                par=c("param1", "param2")
            )
        },
        nbinom={
            Transform(
                loglik_nbinom_base,
                transforms = list(mu~.sensfun2(param1, param2, lambda)),
                par=c("ll.k","param1", "param2")
            )
        },
        nbinom1={
            Transform(
                loglik_nbinom1_base,
                transforms = list(mu~.sensfun2(param1, param2, lambda)),
                par=c("ll.phi","param1", "param2")
            )
        }
    )
    
    hessian_model <- Transform(
        model,
        transforms = list(param1~.trans(param1, transfun_p1, invfun_p1, invfun2_p1), 
                          param2~.trans(param2, transfun_p2, invfun_p2, invfun2_p2)),
        par=model@par,
        keep_hessian=TRUE
    )
    
    hessian_model
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
    
}

##' find Hessian
##' @param data data frame with tvec/count
##' @param params parameter vector
findHess <- function(data, params, 
                     dist=c("gaussian", "poisson", "quasipoisson", "nbinom", "nbinom1")){
    dist <- match.arg(dist)
    t <- data$tvec
    tpars <- trans.pars(params)
    r <- SIR.detsim.hessian(t, tpars)
    
    with(as.list(c(tpars, r)),{
        x <- data$count
        parVec <- c(beta, gamma, N, i0)
        
        sensVec <- c(beta, gamma, N, i0^2*exp(-qlogis(i0)))
        
        sensVec2 <- c(1, 1, 1, 2*i0*exp(-qlogis(i0)) - i0^2 * exp(-qlogis(i0)) * 1/(i0-i0^2))
        
        derVec <- list(nu_I_b, nu_I_g, nu_I_N, nu_I_i)
        
        derMat <- matrix(list(nu_I_bb, nu_I_bg, nu_I_bN, nu_I_bi,
                           nu_I_bg, nu_I_gg, nu_I_gN, nu_I_gi,
                           nu_I_bN, nu_I_gN, nu_I_NN, nu_I_Ni,
                           nu_I_bi, nu_I_gi, nu_I_Ni, nu_I_ii), 4, 4)
        
        findDeriv <- function(i1, i2, dist){
            d1 <- derVec[[i1]]
            d2 <- derVec[[i2]]
            db <- derMat[i1,i2][[1]]
            
            if(i1 == i2){
                d12 <- sensVec2[i1]
            }else{
                d12 <- 0
            }
            ## I don't know if this is right.. I'll have to go through the derivation a few more times...
            if(dist == "norm"){
                n <- length(x)
                sigma2 <- sum((I-x)^2)/(n-1)
                term1 <- sum((I-x)/sigma2 * d1 + (1/(2*sigma2) - ((I - x)^2)/(2*sigma2^2)) * sum(2 * (I - x)/(n-1) * d1)) * d12
                term2 <- (1/sigma2 * d2 - (I-x)/(sigma2^2) * sum(2 * (I - x)/(n-1) * d2)) * d1 + (I - x)/sigma2 * db +
                    (-1/(2 * sigma2^2) * sum(2 * (I - x)/(n-1) * d2) - 
                         (I - x)/(sigma2^2) * d2 + (I - x)^2/(sigma2^3) * sum(2 * (I - x)/(n-1) * d2)) *
                    sum(2 * (I - x)/(n-1) * d1) +
                    (1/(2 * sigma2) - ((I - x)^2)/(2 * sigma2^2)) * sum(2/(n-1) * d2 * d1 + 2 * (I - x)/(n-1) * db)
                deriv <- (term1 + sum(term2) * sensVec[i1]) * sensVec[i2]
            }else if(dist == "pois"){
                deriv <- sum((db + x/I^2 * d1 * d2 - x/I * db) * sensVec[i1] + (1 - x/I) * d1 * d12) * sensVec[i2]
            }

            return(deriv)
        }
        
        m <- matrix(NA, 4, 4)
        
        for(i in 1:4){
            for(j in 1:4){
                m[i,j] <- findDeriv(i,j,dist)
            }
        }
        
        return(m)
    })
    
}
