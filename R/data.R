##' Bombay data set
##'
##' Bombay data: gotten from somewhere (?)
##' Make sure to read the reference for cautions!
##' n.b. fitting obviously isn't working right yet.
##' 
##' @docType data
##' @keywords datasets
##' @name bombay
##' @usage data(bombay)
##' @format A data frame with two elements
##' @references Bacaër, Nicolas. 2012. “The Model of Kermack and McKendrick for the Plague Epidemic in Bombay and the Type Reproduction Number with Seasonality.” Journal of Mathematical Biology 64 (3): 403–22. doi:10.1007/s00285-011-0417-5.
##' @examples
##' par(las=1,bty="l")
##' plot(mort~week,data=bombay)
##' \dontrun{
##' ff <- fitsir(setNames(bombay,c("tvec","count")))
##' ss <- with(bombay,SIR.detsim(t=week,params=trans.pars(coef(ff))))
##' lines(bombay$week,ss)
##' }
NULL
