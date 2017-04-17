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
##' @slot call Objecto of class "\code{language}". 
##' The call that generated the "fitsir" object.
##' @slot coef Object of class "\code{matrix}". 
##' Estimated coefficients and standard errors.
##' @slot summary Object of class "\code{matrix}".
##' Summary of estimated coefficients and standard errors.
##' @slot m2logL Object of class "\code{numeric}".
##' Minus twice the log likelihood.
##' @seealso \code{\link{summary.mle2-class}}
##' @exportClass summary.fitsir
setClass("summary.fitsir", contains="summary.mle2",
         slots= c(call="language",
                  coef="matrix",
                  summary="matrix",
                  m2logL="numeric"))
