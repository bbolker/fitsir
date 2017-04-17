stopifnot(require("testthat"), require("fitsir"))

context("basic tests")
test_that("bombay fits", {
    bombay2 <- setNames(bombay,c("times","count"))
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

test_that("philadelphia fits", {
    phila1918a <- with(phila1918, data.frame(times=seq_along(date), count=pim))
    ss0 <- startfun(data=phila1918a, type="death")
    oldval <- structure(c(1.15432056672458, 0.200572186303702, 0.769401628254234, 
                          6.04738321959361e-06, 0.270409668445559, 44714.8830089917, 44715.1534186602
                          ), .Names = c("R0", "r", "infper", "i0", "I0", "S0", "N"))
    
    expect_equal(summarize.pars(ss0),
                 oldval,
                 tolerance=1e-6)
    
    suppressWarnings(f1 <- fitsir(phila1918a,start=ss0,type="death", method="BFGS"))
    oldval <- structure(c(2.30212921317221, 0.316802420552837, 4.11022494998596, 
                          3.51779353016634e-06, 0.0529188458970781, 15043.1397652272, 15043.1926840731
                          ), .Names = c("R0", "r", "infper", "i0", "I0", "S0", "N"))
    
    expect_equal(summarize.pars(coef(f1)),
                 oldval,
                 tolerance=1e-6)
        
})

test_that("harbin fits", {
    ss0 <- startfun(data=harbin, tcol="week", icol="Deaths", type="death")
    oldval <- structure(c(1.9397042255199, 0.781494541260066, 1.20244502796493, 
                         0.000107170076053567, 0.158254392786069, 1476.50760807211, 1476.66586246489
                         ), .Names = c("R0", "r", "infper", "i0", "I0", "S0", "N"))
    
    expect_equal(summarize.pars(ss0),
                 oldval,
                 tolerance=1e-6)
    
    suppressWarnings(f1 <- fitsir(harbin,start=ss0,type="death",method="BFGS",dist="nbinom1",tcol="week",icol="Deaths"))
    oldval <- structure(c(1.86307149474402, 0.79799269531904, 1.08155312674755, 
                          0.000351206271779594, 0.700632329574636, 1994.23051176552, 1994.93114409509
                          ), .Names = c("R0", "r", "infper", "i0", "I0", "S0", "N"))
    
    expect_equal(summarize.pars(coef(f1)),
                 oldval,
                 tolerance=1e-6)
    
})
