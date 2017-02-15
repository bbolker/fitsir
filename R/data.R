##' Bombay data set
##'
##' Bombay data: gotten from somewhere (?)
##' Make sure to read the reference for cautions!
##' n.b. fitting obviously isn't working right yet.
##' 
##' @docType data
##' @encoding UTF-8
##' @keywords datasets
##' @name bombay
##' @usage data(bombay)
##' @format A data frame with two elements
##' @references Bacaër, Nicolas. 2012. “The Model of Kermack and McKendrick for the Plague Epidemic in Bombay and the Type Reproduction Number with Seasonality.” Journal of Mathematical Biology 64 (3): 403–22. doi:10.1007/s00285-011-0417-5.
##' @examples
##' par(las=1,bty="l")
##' plot(mort~week,data=bombay)
##' \dontrun{
##' ## NONE OF THESE ACTUALLY WORK RELIABLY!
##' ff <- fitsir(setNames(bombay,c("tvec","count")))
##' ff3 <- fitsir(setNames(bombay,c("tvec","count")),start=startfun(log.N=10))
##' ff2 <- fitsir(setNames(bombay,c("tvec","count")),method="BFGS")
##' ff4 <- fitsir(setNames(bombay,c("tvec","count")),start=startfun(log.N=10),
##'         method="SANN")
##' ss <- with(bombay,SIR.detsim(t=week,params=trans.pars(coef(ff3))))

##' lines(bombay$week,ss)
##' }
NULL

##' 1918 Philadelphia flu data set
##' 
##' @docType data
##' @encoding UTF-8
##' @keywords datasets
##' @name phila1918
##' @usage data(phila1918)
##' @format A data frame with two elements
NULL

##' 1912 Harbin plague data set
##' 
##' @docType data
##' @encoding UTF-8
##' @keywords datasets
##' @name harbin
##' @usage data(harbin)
##' @format A data frame with two elements
NULL