
##' gradient function (solves with {S,log(I)} for stability)
##' @param t time vector
##' @param x state vector (with names \code{S}, \code{logI})
##' @param params parameter vector (with names \code{beta} (scaled transmission rate), \code{N} (population size), \code{gamma} (recovery/removal rate)
##' @return gradient vector for a simple SIR model
##' @export
SIR.grad <- function(t, x, params) {
    g <- with(as.list(c(x,params)), {
        I <- exp(logI)
        dS <- -beta*exp(logI)*S/N
        dlogI <- beta*S/N-gamma
        list(c(dS,dlogI), I = I)
    })
}

##' @rdname SIR.grad
SIR.grad.sens <- function(t, x, params) {
    g <- with(as.list(c(x,params)), {
        I <- exp(logI)
        dS <- -beta*S*I/N
        dlogI <- beta*S/N-gamma
                  
        grad_SS <- - beta * I/N
        grad_SI <- - beta * S/N
        grad_IS <- beta*I/N
        grad_II <- beta*S/N-gamma
                  
        dnu_S_b <- grad_SS * nu_S_b + grad_SI * nu_I_b - S*I/N
                  
        dnu_S_N <- grad_SS * nu_S_N + grad_SI * nu_I_N + beta*S*I/N^2
                  
        dnu_S_g <- grad_SS * nu_S_g + grad_SI * nu_I_g
                  
        dnu_S_i <- grad_SS * nu_S_i + grad_SI * nu_I_i
                  
        dnu_I_b <- grad_IS * nu_S_b + grad_II * nu_I_b + S*I/N
                  
        dnu_I_N <- grad_IS * nu_S_N + grad_II * nu_I_N - beta*S*I/N^2
                  
        dnu_I_g <- grad_IS * nu_S_g +  grad_II * nu_I_g - I
                  
        dnu_I_i <- grad_IS * nu_S_i + grad_II * nu_I_i
                  
        list(c(dS, dlogI, dnu_S_b, dnu_S_g, dnu_S_N, dnu_S_i, dnu_I_b, dnu_I_g, dnu_I_N, dnu_I_i), I = I)
    })
}