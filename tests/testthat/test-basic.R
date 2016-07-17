stopifnot(require("testthat"), require("fitsir"))

context("basic tests")
test_that("bombay fits", {
    data("bombay")
    bombay2 <- setNames(bombay,c("tvec","count"))
    ss0 <- startfun(data=bombay2,auto=TRUE)
    expect_equal(summarize.pars(ss0),
                 structure(c(2.65618522, 0.34685520, 4.77486057, 
                             0.00368742173477816, 5.02922941, 1363.88777),
                           .Names = c("R0", "r", "infper", "i0", "I0", "N")),
                 tolerance=1e-6)
    f1 <- fitsir(bombay2,start=ss0)
    expect_equal(summarize.pars(bbmle::coef(f1)),
                 structure(c(1.1951349233, 0.394289728, 0.494902375,
                             6.47960233e-05, 4.013979769, 61947.9338), 
                           .Names = c("R0", "r", "infper", "i0", "I0", "N")),
                 tolerance=1e-6)
})
    
