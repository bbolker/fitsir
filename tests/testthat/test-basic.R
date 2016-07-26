stopifnot(require("testthat"), require("fitsir"))

context("basic tests")
test_that("bombay fits", {
    data("bombay")
    bombay2 <- setNames(bombay,c("tvec","count"))
    ss0 <- startfun(data=bombay2,auto=TRUE)
    expect_equal(summarize.pars(ss0),
                 structure(c(2.91753342360012, 0.346814013295034, 5.5289963787273, 
                             0.00388686417468641, 5.02922941484824, 1288.874850795, 1293.90408020985
                             ), .Names = c("R0", "r", "infper", "i0", "I0", "S0", "N")),

                 tolerance=1e-6)
    f1 <- fitsir(bombay2,start=ss0)
    expect_equal(summarize.pars(coef(f1)),
                 structure(c(1.01448544704402, 0.390081804833677, 0.0371343827487676, 
4.44877097202117e-07, 3.7043246924268, 8326621.14492856, 8326624.84925325
), .Names = c("R0", "r", "infper", "i0", "I0", "S0", "N")),
                 tolerance=1e-6)
})
    
