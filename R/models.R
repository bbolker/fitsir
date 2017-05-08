setClass(
    "loglik.fitsir",
    slots = c(
        name = "character",
        expr = "expression",
        count = "character",
        mean = "character",
        par = "character",
        grad = "list",
        hessian = "list"
    )
)

setMethod(
    "initialize",
    "loglik.fitsir",
    definition = function(.Object, name, model, count="X", mean, par=NULL,
                          keep_grad=TRUE,
                          keep_hessian=FALSE) {
        .Object@name <- name
        if (!is(model, "formula"))
            stop("model must be a formula")
        f <- as.list(model)
        .Object@expr <- as.expression(f[[3]])
        .Object@count <- count
        .Object@mean <- mean
        .Object@par <- par <- as.character(par)
        # compute the gradient
        vars <- c(mean, par)
        deriv <- function(expr) {
            d <- lapply(vars, function(p){Deriv::Deriv(expr, p)})
            names(d) <- vars
            d
        }
        if(keep_grad) {
            .Object@grad <- deriv(.Object@expr)
        } else .Object@grad <- list()
        
        if(keep_hessian) {
            .Object@hessian <- lapply(.Object@grad, function(d) deriv(d))
            names(.Object@hessian) <- vars
        } else .Object@hessian <- list()
        
        
        .Object
    }
)

##' @docType methods
##' @export
setGeneric(
    "Eval",
    def = function(object, ...) {
        standardGeneric("Eval")
    }
)

##' Evaluate a model
##' @param data a dataframe object holding inut values, if NULL, take from ...
##' @param par a named vector (or list) of parameter values, if NULL, take from ...
##' @param ... the input values and parameter values
##' @return numeric
##' @docType methods
##' @export
setMethod(
    "Eval",
    "loglik.fitsir",
    definition = function(object, count, mean, par=NULL, ...) {
        frame <- list(count, mean, par)
        names(frame) <- c(object@count, object@mean, object@par)
        frame <- append(frame, list(...))
        eval(object@expr, frame)
    }
)

##' @docType methods
##' @export
setGeneric(
    "grad",
    def = function(object, ...) {
        standardGeneric("grad")
    }
)

##' the gradients w.r.t. to parameters
##' compute the gradient of the model as a function of the parameters
##' @param data a dataframe object holding inut values, if NULL, take from ...
##' @param par a named vector (or list) of parameters to compute the derivatives
##' @param ... the input values and parameter values
##' @return a dataframe with each column as a partial derivative values
##' @docType methods
##' @export
setMethod(
    "grad",
    "loglik.fitsir",
    definition <- function(object, count, mean, par, ...) {
        frame <- list(count, mean, par)
        names(frame) <- c(object@count, object@mean, object@par[!grepl("param", object@par)])
        frame <- append(frame, list(...))
        var <- c(object@mean, object@par)
        l <- lapply(object@grad[var], function(deriv) { eval(deriv, frame)})
        l
    }
)

#' @docType methods
#' @export
setGeneric(
    "hessian",
    def = function(object, ...) {
        standardGeneric("hessian")
    }
)

##' the hessian w.r.t. to parameters
##' compute the hessian of the model as a function of the parameters
##' @param data a dataframe object holding inut values, if NULL, take from ...
##' @param par a named vector (or list) of parameters to compute the hessians
##' @param ... the input values and parameter values
##' @return a dataframe with each column as a partial derivative values
##' @docType methods
##' @export
setMethod(
    "hessian",
    "loglik.fitsir",
    definition <- function(object, data=NULL, par, ...) {
        frame <- list(count, mean, par)
        names(frame) <- c(object@count, object@mean, object@par)
        var <- c(object@mean, object@par)
        l <- lapply(object@hessian[var], function(dd) { 
            sapply(dd[var], function(e) eval(e, frame))})
        n <- length(l)
        mat <- array(dim=c(dim(l[[1]]), n), dimnames=list(NULL, var, var))
        for (i in 1:n)
            mat[,,i] <- l[[i]]
        mat
    }
)

#' @docType methods
#' @export
setGeneric(
    "Transform",
    def = function(object, ...) {
        standardGeneric("Transform")
    }
)

setMethod(
    "Transform",
    "loglik.fitsir",
    definition <- function(object, transforms=NULL, count="X", mean, par) {
        trans <- function(formulae, allvars) {
            # extract vars from an expression
            vars <- function(e) {
                if (is.numeric(e)) return(c())
                if (is.name(e)) return(as.character(e))
                if (!is.call(e)) 
                    stop("unknown class: ", class(e))
                v <- c()
                for (i in 2:length(e)) {
                    v <- c(v, vars(e[[i]]))
                }
                v
            }
            
            l <- list()
            for (f in formulae) {
                if (!is(f, "formula"))
                    stop("transforms must be formula: ", as.character(f))
                var <- as.character(f[[2]])
                if (!var %in% allvars) next
                input <- vars(f[[3]])
                l[[var]] <- f[[3]]
            }
            l
        }
        # substitute the transformation expressions
        subst <- function(e) {
            if (is.numeric(e)) return(e)
            if (is.name(e)) {
                v <- as.character(e)
                expr <- transforms[[v]]
                if (is.null(expr)) return(e)
                return (expr)
            }
            if (!is.call(e)) stop("unknow class: ", class(e))
            l <- list(e[[1]])
            for (i in 2:length(e))
                l <- c(l, subst(e[[i]]))
            as.call(l)
        }
        # if no transform, return model
        if (length(transforms) == 0) 
            return(object)
        allvars <- c(object@count, object@mean, object@par)
        transforms <- trans(transforms, allvars)
        f <- c(as.symbol("~"), as.symbol("LL"), subst(object@expr[[1]]))
        f <- as.formula(as.call(f))
        
        if (missing(mean)) mean <- object@mean
        if (missing(par)) par <- object@par
        
        new("loglik.fitsir", object@name, f, count, mean=mean, par=par)
    }
)

## Trying to get partial derivatives to work...
.sensfun <- function(x, mean) mean

##' @importFrom Deriv drule
drule[[".sensfun"]] <- alist(x=nu, mean=1)

loglik_gaussian <- new("loglik.fitsir", "gaussian",
                       LL ~ -(X-mu)^2/(2*sigma^2) - log(sigma) - 1/2*log(2*pi),
                       mean="mu", par="sigma")
loglik_gaussian <- Transform(loglik_gaussian,
    transforms = list(sigma ~ sqrt(sum((X-mu)^2)/(length(X)-1)))
)

loglik_gaussian_s <- Transform(
    loglik_gaussian,
    transforms = list(mu~.sensfun(param, mu)), 
    par="param"
)

loglik_poisson <- new("loglik.fitsir", "poisson", 
                      LL ~ X*log(lambda) - lambda - lgamma(X+1), 
                      mean = "lambda", par = c())
loglik_poisson_s <- Transform(
    loglik_poisson,
    transforms = list(lambda~.sensfun(param, lambda)),
    par="param"
)

loglik_nbinom <- new ("loglik.fitsir", "nbinom",
                      LL ~ lgamma(ll.k+X) - lgamma(ll.k) - lgamma(X+1) +
                          ll.k*log(ll.k) - ll.k*log(ll.k+mu) + X*log(mu) - 
                          X*log(ll.k+mu),
                      mean="mu",
                      par = "ll.k")
loglik_nbinom <- Transform(
    loglik_nbinom,
    transforms = list(ll.k ~ exp(ll.k))
)
loglik_nbinom_s <- Transform(
    loglik_nbinom,
    transforms = list(mu~.sensfun(param, mu)),
    par=c("ll.k", "param")
)

# negative binomial '1' likelihood
# var proportional to mean
# v=mu*(1+mu/k), k>0
# v=mu*(1+phi), phi>0
# i.e. mu/k=phi -> k=mu/phi
#' @export
loglik_nbinom1 <- new ("loglik.fitsir", "nbinom",
                       LL ~ lgamma(mu/ll.phi+X) - lgamma(mu/ll.phi) - lgamma(X+1) +
                           mu/ll.phi*log(mu/ll.phi) - mu/ll.phi*log(mu/ll.phi+mu) + X*log(mu) - 
                           X*log(mu/ll.phi+mu),
                       mean = "mu",
                       par = "ll.phi")
loglik_nbinom1 <- Transform(
    loglik_nbinom1,
    transforms = list(ll.phi ~ exp(ll.phi))
)
loglik_nbinom1_s <- Transform(
    loglik_nbinom1,
    transforms = list(mu~.sensfun(param, mu)),
    par=c("ll.phi", "param")
)
