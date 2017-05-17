fitsirx <- function(data, start,
                     dist,
                     type) {
    
    f.env <- new.env()
    ## set initial values
    assign("oldnll",NULL,f.env)
    assign("oldpar",NULL,f.env)
    assign("oldgrad",NULL,f.env)
    assign("data", data, f.env)
    objfun <- function(par, count, times, dist, type) {
        if (identical(par,oldpar)) {
            return(oldnll)
        }
        
        v <- SIR.sensitivity(par, count, times, dist, type)
        oldnll <<- v[1]
        oldgrad <<- v[-1]
        oldpar <<- par
        
        return(oldnll)
    }
    environment(objfun) <- f.env
    gradfun <- function(par, count, times, dist, type) {
        if (identical(par,oldpar)) {
            return(oldgrad)
        }
        v <- SIR.sensitivity(par, count, times, dist, type)
        oldnll <<- v[1]
        oldgrad <<- v[-1]
        oldpar <<- par
        return(oldgrad)
    }
    environment(gradfun) <- f.env
    fit <- optimx(start, objfun, gr = gradfun,
                  control = list(all.methods = TRUE, trace = 0),
                  count=data$count,
                  times=data$times,
                  type=type,
                  dist=dist)
    return(fit)
}