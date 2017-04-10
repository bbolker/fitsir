stopifnot(require("testthat"), require("fitsir"))

context("basic tests")
test_that("bombay fits", {
    data("bombay")
    bombay2 <- setNames(bombay,c("tvec","count"))
    ss0 <- startfun(data=bombay2)
    oldval <- structure(c(1.0309058069171, 0.346654710263127, 0.0891544410103143, 
                          2.2990741991925e-06, 4.0426380952539, 1758372.47982203, 1758376.52246012
                          ), .Names = c("R0", "r", "infper", "i0", "I0", "S0", "N"))

    expect_equal(summarize.pars(ss0),
                 oldval,
                 tolerance=1e-6)
    f1 <- fitsir(bombay2,start=ss0)
    oldval <- structure(c(1.03291670446211, 0.379613660510302, 0.086711064132579, 
2.72670847076781e-06, 4.58863236589658, 1682842.11650311, 1682846.70513547
), .Names = c("R0", "r", "infper", "i0", "I0", "S0", "N"))

    expect_equal(summarize.pars(coef(f1)),
              oldval,
                 tolerance=1e-6)
})
    
