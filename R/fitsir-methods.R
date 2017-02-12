setMethod("plot", signature(x = "fitsir.mle2", y = "missing"),
    function(x, ...){
        with(as.list(x@data), {
            pars <- coef(x)
            i.hat <- SIR.detsim(tvec, trans.pars(pars), type)
            plot(tvec, i.hat, type = "l", ...)
            points(tvec, count)
            
            invisible()
        })
    }
)