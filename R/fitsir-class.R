##' Class "fitsir". 
##' Result of SIR model fitting based on Maximum Likelihood Estimation
##' 
##' @name fitsir-class
##' @rdname fitsir-class
##' @seealso \code{\link{mle2-class}}
##' @exportClass fitsir
setClass("fitsir", contains="mle2")

setClass("summary-fitsir", representation(call="language",
                                          coef="vector",
                                          summary="matrix",
                                          m2logL="numeric"))