## let's say we have findSens defined
findSens <- function(p) {
    return(c(SSQ=sum(p^2),p^2))
}
f.env <- new.env()
objfun <- function(par) {
   if (identical(par,oldpar)) {
      return(oldSSQ)
   }
   v <- findSens(par)
   oldSSQ <<- v["SSQ"]
   oldgrad <<- v[-1]
   oldpar <<- par
   return(oldSSQ)
}
environment(objfun) <- f.env
gradfun <- function(par) {
   if (identical(par,oldpar)) {
      return(oldsens)
   }
   v <- findSens(par)
   oldSSQ <<- v["SSQ"]
   oldgrad <<- v[-1]
   oldpar <<- par
   return(oldsens)
}
environment(gradfun) <- f.env

a <- 5
f <- function() {
     r <- a^2
     a <- 1
     a <<- 10
     r <- a^2
     return(r)
}

```
