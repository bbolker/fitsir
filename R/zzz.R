

.onLoad <- function(libname, pkgname) {
    ## prevent NOTE false positives from R CMD check ...
    N <- a <- b <- i0 <- nu_I_N <- nu_I_NN <- nu_I_Ni <- nu_I_b <-
        nu_I_bN <- nu_I_bb <- nu_I_bg <- NULL

    drule[["lbeta"]] <- drule[["w_lbeta"]] <- alist(a=dfun(a,b),
                                                    b=dfun(b,a))
    drule[["dfun"]] <- alist(x=dfun2(x,y),
                             y=dfun2(y,x))
    
    drule[[".sensfun"]] <- alist(beta=nu_I_b, gamma=nu_I_g, N=nu_I_N, i0=nu_I_i, mean=1)
    
    drule[[".sensfun2"]] <- alist(beta=.nu_beta(beta,gamma,N,i0, nu_I_b),
                                  gamma=.nu_gamma(beta,gamma,N,i0, nu_I_g),
                                  N=.nu_N(beta,gamma,N,i0, nu_I_N),
                                  i0=.nu_i(beta,gamma,N,i0, nu_I_i))
    
    drule[[".nu_beta"]] <- alist(beta=nu_I_bb, gamma=nu_I_bg, N=nu_I_bN, i0=nu_I_bi)
    
    drule[[".nu_gamma"]] <- alist(beta=nu_I_bg, gamma=nu_I_gg, N=nu_I_gN, i0=nu_I_gi)
    
    drule[[".nu_N"]] <- alist(beta=nu_I_bN, gamma=nu_I_gN, N=nu_I_NN, i0=nu_I_Ni)
    
    drule[[".nu_i"]] <- alist(beta=nu_I_bi, gamma=nu_I_gi, N=nu_I_Ni, i0=nu_I_ii)
}
