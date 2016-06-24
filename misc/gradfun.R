## let's say we have findSens defined
## (this is silly but works)
findSens <- function(p) {
    return(c(SSQ=sum(p^2),p^2))
}
f.env <- new.env()
## set initial values
assign("oldpar",NULL,f.env)
assign("oldgrad",NULL,f.env)
objfun <- function(par,verbose=TRUE) {
   if (identical(par,oldpar)) {
       if (verbose) cat("returning old version of SSQ\n")
       return(oldSSQ)
   }
   if (verbose) cat("computing new version (SSQ)\n")
   v <- findSens(par)
   oldSSQ <<- v["SSQ"]
   oldgrad <<- v[-1]
   oldpar <<- par
   return(oldSSQ)
}
environment(objfun) <- f.env
gradfun <- function(par,verbose=TRUE) {
   if (identical(par,oldpar)) {
       if (verbose) cat("returning old version of grad\n")
      return(oldgrad)
   }
   if (verbose) cat("computing new version (grad)\n")
   v <- findSens(par)
   oldSSQ <<- v["SSQ"]
   oldgrad <<- v[-1]
   oldpar <<- par
   return(oldgrad)
}
environment(gradfun) <- f.env

## example ...
objfun(1:3)
gradfun(1:3)
objfun(2:4)
gradfun(2:4)
