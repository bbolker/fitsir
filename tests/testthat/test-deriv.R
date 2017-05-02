stopifnot(require("testthat"), require("fitsir"))

context("derivative test")
test_that("sensitivity", {
    ss <- startfun()
    
    SIR.detsim2 <- function(t, params, type) SIR.detsim(t, trans.pars(params), type)
    
    for (type in c("prevalence", "incidence", "death")) {
        expect_equal(
            numDeriv::jacobian(SIR.detsim2, startfun(), t=1:10, type=type),
            unname(as.matrix(SIR.detsim(1:10, trans.pars(startfun()), grad=TRUE, type=type)[,-1])),
            tolerance=1e-5
        )
    }
})

test_that("likelihood", {
    ss <- startfun()
    
    for (type in c("prevalence", "incidence", "death")) {
        for(dist in c("gaussian", "poisson")) {
            expect_equal(
                numDeriv::grad(SIR.logLik, ss, count=harbin$Deaths, times=harbin$week, type=type, dist=dist),
                unname(SIR.sensitivity(ss, harbin$Deaths, harbin$week, type=type, dist=dist)[-1])
            )
        }
    }
    
    ss <- c(startfun(), 5)
    
    for (type in c("prevalence", "incidence", "death")) {
        for(dist in c("nbinom", "nbinom1")) {
            expect_equal(
                numDeriv::grad(SIR.logLik, ss, count=harbin$Deaths, times=harbin$week, type=type, dist=dist),
                unname(SIR.sensitivity(ss, harbin$Deaths, harbin$week, type=type, dist=dist)[-1])
            )
        }
    }
})
