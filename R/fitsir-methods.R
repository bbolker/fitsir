setMethod("plot", signature(x = "fitsir.mle2", y = "missing"),
    function(x, type = "l", ...){
        with(as.list(x@data), {
            pars <- coef(x)
            i.hat <- SIR.detsim(tvec, trans.pars(pars), incidence = incidence)
            plot(tvec, i.hat, type = type, ...)
            points(tvec, count)
            
            invisible()
        })
    }
)