stopifnot(require("testthat"), require("fitsir"))

context("derivative test")
test_that("sensitivity", {
    ss <- trans.pars(startfun(harbin, tcol="week", icol="Deaths"))
    
    for (type in c("prevalence", "incidence", "death")) {
        expect_equal(
            numDeriv::jacobian(SIR.detsim, ss, t=1:10, type=type),
            unname(as.matrix(SIR.detsim(1:10, ss, grad=TRUE, type=type)[,-1])),
            tolerance=1e-5
        )
    }
})
