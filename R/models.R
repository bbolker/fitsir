##' Class representing log-likelihood models used to fit the SIR model
##' 
##' @slot name name of the distribution
##' @slot expr an expressioin specifying the model
##' @slot count observation variable name
##' @slot mean mean variable name
##' @slot par additional parameter names
##' @slot grad the gradient with respect to the parameters
##' @slot hessian the Hessian of the model with respect 
##' @exportClass loglik.fitsir
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

##' the initializer for loglik.fitsir
##' 
##' @param name name of the distribution
##' @param model the formula specifying the model
##' @param count observation variable name
##' @param mean mean variable name
##' @param par additional parameter names
##' @param keep_grad maintain the gradient as part of the model
##' @param keep_hessian maintain the hessian as part of the model
##' @docType methods
##' @exportMethod initialize
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
        vars <- c(par)
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

setGeneric(
    "Eval",
    def = function(object, ...) {
        standardGeneric("Eval")
    }
)

##' Evaluate a model
##' @param object object to be evaluated
##' @param count observations
##' @param mean mean values
##' @param par additional parameters
##' @param ... other values if required
##' @return numeric
##' @docType methods
##' @exportMethod Eval
setMethod(
    "Eval",
    "loglik.fitsir",
    definition = function(object, count, mean, par=NULL, ...) {
        frame <- list(count, mean, par)
        names(frame) <- c(object@count, object@mean, object@par[!grepl("param", object@par)])
        frame <- append(frame, list(...))
        eval(object@expr, frame)
    }
)

setGeneric(
    "grad",
    def = function(object, ...) {
        standardGeneric("grad")
    }
)

##' Evaluate the gradient of a model
##' @param object object to be evaluated
##' @param count observations
##' @param mean mean values
##' @param par additional parameters
##' @param var character vector specifying names of partial derivative parameters
##' @param ... other values if required
##' @return a list with each element as a partial derivative values
##' @docType methods
##' @exportMethod grad
setMethod(
    "grad",
    "loglik.fitsir",
    definition <- function(object, count, mean, par, var, ...) {
        frame <- list(count, mean)
        frame <- append(frame, par)
        names(frame) <- c(object@count, object@mean, object@par[!grepl("param", object@par)])
        frame <- append(frame, list(...))
        if (missing(var)) var <- c(object@par)
        l <- lapply(object@grad[var], function(deriv) { eval(deriv, frame)})
        l
    }
)

setGeneric(
    "hessian",
    def = function(object, ...) {
        standardGeneric("hessian")
    }
)

##' Evaluate the hessian of a model
##' @param object object to be evaluated
##' @param count observations
##' @param mean mean values
##' @param par additional parameters
##' @param var character vector specifying names of partial derivative parameters
##' @param ... other values if required
##' @return a list with each element as a partial derivative values
##' @docType methods
##' @exportMethod hessian
setMethod(
    "hessian",
    "loglik.fitsir",
    definition <- function(object, count, mean, par, var, ...) {
        frame <- list(count, mean)
        frame <- append(frame, par)
        names(frame) <- c(object@count, object@mean, object@par[!grepl("param", object@par)])
        frame <- append(frame, list(...))
        if (missing(var)) var <- c(object@par)
        l <- lapply(object@hessian[var], function(dd) { 
            sapply(dd[var], function(e) eval(e, frame))})
        n <- length(l)
        mat <- array(dim=c(dim(l[[1]]), n), dimnames=list(NULL, var, var))
        for (i in 1:n)
            mat[,,i] <- l[[i]]
        mat
    }
)

setGeneric(
    "Transform",
    def = function(object, ...) {
        standardGeneric("Transform")
    }
)

##' Transform the model
##' @param object object
##' @param transforms list of formulas specifying transformations
##' @param count observation variable name
##' @param mean mean variable name
##' @param par additional parameter names
##' @param keep_grad maintain the gradient as part of the model
##' @param keep_hessian maintain the hessian as part of the model
##' @return loglik.fitsir object
##' @docType methods
##' @exportMethod Transform
setMethod(
    "Transform",
    "loglik.fitsir",
    definition <- function(object, transforms=NULL, 
                           name,
                           count="X", 
                           mean, par,
                           keep_grad=TRUE,
                           keep_hessian=FALSE) {
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
        
        if (missing(name)) name <- object@name
        if (missing(mean)) mean <- object@mean
        if (missing(par)) par <- object@par
        
        new("loglik.fitsir", name, f, count, mean=mean, par=par, keep_grad=keep_grad, keep_hessian=keep_hessian)
    }
)

## Trying to get partial derivatives to work...
.sensfun <- function(x, mean) mean

.sensfun2 <- function(beta, gamma, N, i0, mean) mean
.nu_beta <- function(beta, gamma, N, i0, nu_I_b) nu_I_b
.nu_gamma <- function(beta, gamma, N, i0, nu_I_g) nu_I_g
.nu_N <- function(beta, gamma, N, i0, nu_I_N) nu_I_N
.nu_i <- function(beta, gamma, N, i0, nu_I_i) nu_I_i

## use Taylor expansion of digamma(a+b) for a>>b
## discontinuity in second derivative, but ... probably OK
##' @importFrom Deriv drule
dfun <- function(x,y,mag=1e8) {
    return(ifelse(x/y>mag,
                  -y*trigamma(x),
                  digamma(x)-digamma(x+y)))
}
dfun2 <- function(x,y,mag=1e8,focal="x") {
    return(switch(focal,
                  x=ifelse(x/y>mag,
                           -y*psigamma(x,2),
                           trigamma(x)-trigamma(x+y)),
                  y=ifelse(x/y>mag,
                           -trigamma(x),
                           -trigamma(x+y))))
}

w_lbeta <- function(a,b) {
    ## when we have an effectively-Poisson casee
    ## lbeta gives "underflow occurred in 'lgammacor'" frequently ...
    ## suppressWarnings() causes an obscure error ?
    ## using w_lbeta rather than lbeta causes obscure errors from Deriv()
    op <- options(warn=-1)
    on.exit(options(op))
    return(lbeta(a,b))
}
drule[["lbeta"]] <- drule[["w_lbeta"]] <- alist(a=dfun(a,b),
                                                b=dfun(b,a))
drule[["dfun"]] <- alist(x=dfun2(x,y),
                         y=dfun2(y,x))

drule[[".sensfun"]] <- alist(x=nu, mean=1)

drule[[".sensfun2"]] <- alist(beta=.nu_beta(beta,gamma,N,i0, nu_I_b),
                              gamma=.nu_gamma(beta,gamma,N,i0, nu_I_g),
                              N=.nu_N(beta,gamma,N,i0, nu_I_N),
                              i0=.nu_i(beta,gamma,N,i0, nu_I_i))

drule[[".nu_beta"]] <- alist(beta=nu_I_bb, gamma=nu_I_bg, N=nu_I_bN, i0=nu_I_bi)

drule[[".nu_gamma"]] <- alist(beta=nu_I_bg, gamma=nu_I_gg, N=nu_I_gN, i0=nu_I_gi)

drule[[".nu_N"]] <- alist(beta=nu_I_bN, gamma=nu_I_gN, N=nu_I_NN, i0=nu_I_Ni)

drule[[".nu_i"]] <- alist(beta=nu_I_bi, gamma=nu_I_gi, N=nu_I_Ni, i0=nu_I_ii)

##' gaussian log-likelihood base model
loglik_gaussian_base<- new("loglik.fitsir", "gaussian",
                       LL ~ -(X-mu)^2/(2*sigma^2) - log(sigma) - 1/2*log(2*pi),
                       mean="mu", par="sigma")

loglik_gaussian_base <- Transform(loglik_gaussian_base,
    transforms = list(sigma ~ sqrt(sum((X-mu)^2)/(length(X)-1)))
)

##' gaussian log-liklihood model with sensitivity
##' @export
loglik_gaussian <- Transform(
    loglik_gaussian_base,
    transforms = list(mu~.sensfun(param, mu)), 
    par="param"
)

##' poisson log-liklihood base model
loglik_poisson_base <- new("loglik.fitsir", "poisson", 
                      LL ~ X*log(lambda) - lambda - lgamma(X+1), 
                      mean = "lambda", par = c())

##' poisson log-liklihood model with sensitivity
##' @export
loglik_poisson <- Transform(
    loglik_poisson_base,
    transforms = list(lambda~.sensfun(param, lambda)),
    par="param"
)

##' negative binomial log-liklihood base model
loglik_nbinom_base <- new ("loglik.fitsir", "nbinom",
                      LL ~ -lbeta(ll.k, X) - log(X) + ll.k * (-log1p(mu/ll.k)) + 
                          X * log(mu) - X * log(ll.k + mu),
                      mean="mu",
                      par = "ll.k")

# negative binomial '1' likelihood
# var proportional to mean
# v=mu*(1+mu/k), k>0
# v=mu*(1+phi), phi>0
# i.e. mu/k=phi -> k=mu/phi
##' negative binomial 1 log-liklihood base model
loglik_nbinom1_base <- Transform(
    loglik_nbinom_base,
    transforms=list(ll.k~mu/ll.phi),
    name="nbinom1",
    par="ll.phi"
)

loglik_nbinom_base <- Transform(
    loglik_nbinom_base,
    transforms = list(ll.k ~ exp(ll.k))
)

##' negative binomial log-liklihood model with sensitivity
##' @export
loglik_nbinom <- Transform(
    loglik_nbinom_base,
    transforms = list(mu~.sensfun(param, mu)),
    par=c("ll.k", "param")
)

loglik_nbinom1_base <- Transform(
    loglik_nbinom1_base,
    transforms = list(ll.phi ~ exp(ll.phi))
)

##' negative binomial 1 log-liklihood model with sensitivity
##' @export
loglik_nbinom1 <- Transform(
    loglik_nbinom1_base,
    transforms = list(mu~.sensfun(param, mu)),
    par=c("ll.phi", "param")
)
