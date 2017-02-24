setMethod("plot", signature(x = "fitsir.mle2", y = "missing"),
    function(x, ...){
        with(as.list(x@data), {
            pars <- coef(x)
            i.hat <- SIR.detsim(tvec, trans.pars(pars), type)
            ymin <- min(i.hat, count)
            ymax <- max(i.hat, count)
            plot(tvec, i.hat, type = "l", ylim = c(ymin, ymax), ...)
            points(tvec, count)
            
            invisible()
        })
    }
)