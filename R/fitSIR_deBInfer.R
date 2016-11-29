SIR.deBInfer <- function(t, x, params){
    g <- with(as.list(c(x,params)),{
        
        gamma <- r/(R0 - 1)
        beta <- r + gamma
        
        inf <- beta*S*I/(S+I+R)
        dS <- -inf
        dI <- inf-gamma*I
        dR <- gamma*I
        list(c(dI,dS,dR))
    })
}

deBInfer.logLik <- function(data, sim.data, samp){
    r <- sum(dnorm2(data$count, sim.data[,"I"], log = TRUE))
    return(r)
}

fitsir.deBInfer <- function(data,
                     R0.args = list(),
                     r.args = list(),
                     S.args = list(),
                     I.args = list(),
                     iter = 10000,
                     plot = FALSE,
                     cnt = 500, ...){
    tvec <- data$tvec
    count <- data$count
    
    default <- list(fixed = FALSE, samp.type = "rw", prior = "unif")
    values <- list(R0=2,r=1,I=10,S=4e5)
    hypers <- list(R0=list(min=1,max=10),
                   r=list(min=0.5,max=5),
                   I=list(min=1,max=20),
                   S=list(min=1e4,max=1e8)) ## hyperparameters
    prop.var <- list(R0=2,r=1,I=1,S=1000)
    cfun <- function(value,hypers,prop.var) {
        append(list(value=value,hypers=hypers,prop.var=prop.var),default)
    }
    defaultlist <- Map(cfun, values, hypers, prop.var)
    
    argslist <- list(R0.args, r.args, I.args, S.args)
    
    for(i in 1:4){
        def <- defaultlist[[i]]
        new <- argslist[[i]]
        
        m <- match(names(new), names(def))
        
        defaultlist[[i]] <- replace(def, m, new)
    }
    
    R <- debinfer_par(name = "R", var.type = "init", fixed = TRUE, value = 0)
    
    dfun <- function(x,name,var.type) {
        do.call(debinfer_par,
                args=append(x,list(name=name,var.type=var.type)))
    }
    debsetupList <- Map(dfun,
                        defaultlist,
                        list("R0", "r", "I", "S"),
                        list("de", "de", "init", "init"))
    
    mcmc.pars <- do.call("setup_debinfer", append(debsetupList, list(R)))
    
    mcmc_samples <- de_mcmc(N = iter, data=data, de.model=SIR.deBInfer,
                            obs.model=deBInfer.logLik, all.params=mcmc.pars,
                            Tmax = max(tvec), data.times=tvec, cnt=cnt,
                            plot=plot, solver="ode", ...)
    class(mcmc_samples) <- append("fitsir_deBInfer", class(mcmc_samples))
    return(mcmc_samples)
}

plot.fitsir_deBInfer <- function(x, plot.type = "trajectory",
                                 burnin = 1000, n = 500,
                                 prob = 0.95,
                                 ...){
    if(plot.type == "trajectory"){
        tvec <- x$data$tvec
        post_traj <- post_sim(x, n = n, times = tvec, burnin = burnin, output = "all", prob = prob)
        plot(x$data, ...)
        lines(tvec, post_traj[["median"]]$I)
        matlines(tvec, post_traj[["HDI"]]$I, col = 2, lty = 2)
    }else{
        deBInfer:::plot.debinfer_result(x, plot.type, ...)
    }
}
