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
findHess <- function(data, params, dist = c("norm", "pois")){
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
