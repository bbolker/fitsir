
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

## second derivs of gradient wrt parameters
SIR.grad.hessian <- function(t, x, params) {
    g <- with(as.list(c(x,params)),
              {
                  I <- exp(logI)
                  dS <- -beta*S*I/N
                  dlogI <- beta*S/N-gamma
                  
                  grad_SS <- - beta * I/N
                  grad_SI <- - beta * S/N
                  grad_IS <- beta*I/N
                  grad_II <- beta*S/N-gamma
                  
                  ##Sens for beta 
                  dnu_S_b <- grad_SS * nu_S_b + grad_SI * nu_I_b - S*I/N
                  ##Need 4 more equations
                  dnu_S_bb <- (- I/N - beta/N * nu_I_b) * nu_S_b + grad_SS * nu_S_bb +  (- S/N - beta/N * nu_S_b) * nu_I_b + grad_SI * nu_I_bb - 
                      I/N * nu_S_b - S/N * nu_I_b
                  
                  dnu_S_bg <- (-beta/N * nu_I_g) * nu_S_b + grad_SS * nu_S_bg + (-beta/N * nu_S_g) * nu_I_b + grad_SI * nu_I_bg -
                      I/N * nu_S_g - S/N * nu_I_g
                  
                  dnu_S_bN <- (beta*I/N^2 - beta/N * nu_I_N) * nu_S_b + grad_SS * nu_S_bN + (beta*S/N^2 - beta/N * nu_S_N) * nu_I_b + grad_SI * nu_I_bN +
                      S * I/N^2 - I/N * nu_S_N - S/N * nu_I_N
                  
                  dnu_S_bi <- (-beta/N * nu_I_i) * nu_S_b + grad_SS * nu_S_bi + (-beta/N * nu_S_i) * nu_I_b + grad_SI * nu_I_bi -
                      I/N * nu_S_i - S/N * nu_I_i
                  
                  ##Sens for gamma
                  dnu_S_g <- grad_SS * nu_S_g + grad_SI * nu_I_g
                  ##Need 3 more equations
                  dnu_S_gg <- (-beta/N * nu_I_g) * nu_S_g + grad_SS * nu_S_gg + (-beta/N * nu_S_g) * nu_I_g + grad_SI * nu_I_gg
                  
                  dnu_S_gN <- (beta*I/N^2 - beta/N * nu_I_N) * nu_S_g + grad_SS * nu_S_gN + (beta*S/N^2 - beta/N * nu_S_N) * nu_I_g + grad_SI * nu_I_gN
                  
                  dnu_S_gi <- (-beta/N * nu_I_i) * nu_S_g + grad_SS * nu_S_gi + (-beta/N * nu_S_i) * nu_I_g + grad_SI * nu_I_gi
                  
                  ##Sens for N
                  dnu_S_N <- grad_SS * nu_S_N + grad_SI * nu_I_N + beta*S*I/N^2
                  ##Need 2 more equations
                  ##dnu_S_NN <- (beta*I/N^2 - beta/N * nu_I_N) * nu_S_N + grad_SS * nu_S_NN + (beta*S/N^2 - beta/N * nu_S_N) * nu_I_N + grad_SI * nu_I_NN -
                  ##    2*beta*S*I/N^3 + beta*S/N^2 * nu_I_N + beta*I/N^2 * nu_S_N
                  dnu_S_NN <- 0
                  
                  dnu_S_Ni <- (-beta/N * nu_I_i) * nu_S_N + grad_SS * nu_S_Ni + (-beta/N * nu_S_i) * nu_I_N + grad_SI * nu_I_Ni +
                      beta*S/N^2 * nu_I_i + beta*I/N^2 * nu_S_i
                  
                  ##Sens for i0
                  dnu_S_i <- grad_SS * nu_S_i + grad_SI * nu_I_i
                  ##Need 1 more equaion
                  dnu_S_ii <- (-beta/N * nu_I_i) * nu_S_i + grad_SS * nu_S_ii + (-beta/N * nu_S_i) * nu_I_i + grad_SI * nu_I_ii
                  
                  ##Sens for beta
                  dnu_I_b <- grad_IS * nu_S_b + grad_II * nu_I_b + S*I/N
                  ##Need 4 more equations
                  dnu_I_bb <- (I/N + beta/N * nu_I_b) * nu_S_b + grad_IS * nu_S_bb + (S/N + beta/N * nu_S_b) * nu_I_b + grad_II * nu_I_bb +
                      I/N * nu_S_b + S/N * nu_I_b
                  
                  dnu_I_bg <- (beta/N * nu_I_g) * nu_S_b + grad_IS * nu_S_bg + (-1 + beta/N * nu_S_g) * nu_I_b + grad_II * nu_I_bg + 
                      I/N * nu_S_g + S/N * nu_I_g
                  
                  dnu_I_bN <- (-beta*I/N^2 + beta/N * nu_I_N) * nu_S_b + grad_IS * nu_S_bN + (-beta*S/N^2 + beta/N * nu_S_N) * nu_I_b + grad_II * nu_I_bN -
                      S * I/N^2 + I/N * nu_S_N + S/N * nu_I_N
                  
                  dnu_I_bi <- (beta/N * nu_I_i) * nu_S_b + grad_IS * nu_S_bi + (beta/N * nu_S_i) * nu_I_b + grad_II * nu_I_bi +
                      I/N * nu_S_i + S/N * nu_I_i
                  
                  ##Sens for gamma
                  dnu_I_g <- grad_IS * nu_S_g + grad_II * nu_I_g - I
                  ##Need 3 more equations
                  dnu_I_gg <- (beta/N * nu_I_g) * nu_S_g + grad_IS * nu_S_gg + (-1 + beta/N * nu_S_g) * nu_I_g + grad_II * nu_I_gg -
                      nu_I_g
                  
                  dnu_I_gN <- (-beta*I/N^2 + beta/N * nu_I_N) * nu_S_g + grad_IS * nu_S_gN + (-beta*S/N^2 + beta/N * nu_S_N) * nu_I_g + grad_II * nu_I_gN -
                      nu_I_N
                  
                  dnu_I_gi <- (beta/N * nu_I_i) * nu_S_g + grad_IS * nu_S_gi + (beta/N * nu_S_i) * nu_I_g + grad_II * nu_I_gi -
                      nu_I_i
                  
                  ##Sens for N
                  dnu_I_N <- grad_IS * nu_S_N + grad_II * nu_I_N - beta*S*I/N^2
                  ##Need 2 more equations
                  ##dnu_I_NN <- (-beta*I/N^2 + beta/N * nu_I_N) * nu_S_N + grad_IS * nu_S_NN + (-beta*S/N^2 + beta/N * nu_S_N) * nu_I_N + grad_II * nu_I_NN +
                  ##    2*beta*S*I/N^3 - beta*S/N^2 * nu_I_N - beta*I/N^2 * nu_S_N
                  dnu_I_NN <- 0
                  
                  dnu_I_Ni <- (beta/N * nu_I_i) * nu_S_N + grad_IS * nu_S_Ni + (beta/N * nu_S_i) * nu_I_N + grad_II * nu_I_Ni -
                      beta*S/N^2 * nu_I_i - beta*I/N^2 * nu_S_i
                  
                  ##Sens for i0
                  dnu_I_i <- grad_IS * nu_S_i + grad_II * nu_I_i
                  ##Need 1 more equaion
                  dnu_I_ii <- (beta/N * nu_I_i) * nu_S_i + grad_IS * nu_S_ii + (beta/N * nu_S_i) * nu_I_i + grad_II * nu_I_ii
                  
                  
                  list(c(dS, dlogI,
                         dnu_S_b, dnu_S_bb, dnu_S_bg, dnu_S_bN, dnu_S_bi,
                         dnu_S_g, dnu_S_gg, dnu_S_gN, dnu_S_gi,
                         dnu_S_N, dnu_S_NN, dnu_S_Ni,
                         dnu_S_i, dnu_S_ii,
                         dnu_I_b, dnu_I_bb, dnu_I_bg, dnu_I_bN, dnu_I_bi,
                         dnu_I_g, dnu_I_gg, dnu_I_gN, dnu_I_gi,
                         dnu_I_N, dnu_I_NN, dnu_I_Ni,
                         dnu_I_i, dnu_I_ii), I = I)
              })
}
