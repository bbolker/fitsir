##' Class "fitsir". 
##' Result of SIR model fitting based on Maximum Likelihood Estimation
##' 
##' @name fitsir-class
##' @rdname fitsir-class
##' @seealso \code{\link{mle2-class}}
##' @exportClass fitsir
setClass("fitsir", contains="mle2")


##' Class "summary.fitsir".
##' Summary of SIR model fit
##' 
##' @name summary.fitsir-class
##' @rdname summary.fitsir-class
##' @method show \code{signature(object="fitsir")}: Show object.
##' @exportClass summary.fitsir
setClass("summary.fitsir", contains="summary.mle2",
         slots= c(call="language",
                  coef="matrix",
                  summary="matrix",
                  m2logL="numeric"))
