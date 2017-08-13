stopifnot(require("testthat"), require("numDeriv"), require("fitsir"))

context("NB test")
loglik_nbinom <- select_model("nbinom")
loglik_nbinom1 <- select_model("nbinom1")

my_dnbinom <- function(X,mu,k) {
    Eval(loglik_nbinom, mean=mu, count=X, par=list(ll.k=log(k)))
}

my_dnbinom1 <- function(X,mu,phi) {
    Eval(loglik_nbinom1, mean=mu, count=X, par=list(ll.phi=log(phi)))
}

dnbinom1 <- function(X,mu,phi,log) dnbinom(x=X, mu=mu, size=mu/phi, log=log)

grad_dnbinom <- function(X,mu,k) {
    fitsir::grad(loglik_nbinom, mean=mu, count=X, par=list(ll.k=log(k)))
}

grad_dnbinom1 <- function(X,mu,phi) {
    fitsir::grad(loglik_nbinom1, mean=mu, count=X, par=list(ll.phi=log(phi)))
}

test_that("NBconst works", {
    expect_equal(my_dnbinom(1, 1, 2), dnbinom(1, mu=1, size=2, log=TRUE))
    expect_equal(my_dnbinom1(1, 1, 2), dnbinom1(1, mu=1, phi=2, log=TRUE))
})

test_that("NBconst works for 0", {
    expect_equal(my_dnbinom(0, 2, 3), dnbinom(0, mu=2, size=3, log=TRUE))
    expect_equal(my_dnbinom1(0, 1, 2), dnbinom1(0, mu=1, phi=2, log=TRUE))
}) 

test_that("NBconst derivative works", {
    expect_equal(
        grad_dnbinom(1, 1, 2)$ll.k,
        numDeriv::grad(function(x) dnbinom(1, mu=1, size=exp(x), log=TRUE), log(2))
    )
    
    expect_equal(
        grad_dnbinom1(1, 1, 2)$ll.phi,
        numDeriv::grad(function(x) dnbinom1(1, mu=1, phi=exp(x), log=TRUE), log(2))
    )
})

test_that("fits containing 0", {
    harbin2 <- setNames(harbin, c("times", "count"))
    harbin2$count[2] <- 0
    ss <- startfun(harbin2, type="death")
    
    ff <- fitsir(harbin2, c(ss, ll.k=2), dist="nbinom")
    oldval <- structure(c(0.4548345079452, -0.388927930093911, 7.07817670051662, 
                          -7.69379344163717, 3.65629428597238), 
                        .Names = c("log.beta", "log.gamma", "log.N", "logit.i", "ll.k"))
    
    expect_equal(
        coef(ff),
        oldval
    )
})
