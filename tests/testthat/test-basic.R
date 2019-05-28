stopifnot(require("testthat"), require("fitsir"))

context("basic tests")
test_that("bombay fits", {
    bombay2 <- setNames(bombay,c("times","count"))
    ss0 <- startfun(data=bombay2)
    oldval <- structure(c(1.02925789958946, 0.346654710263133, 0.0844006982257594, 
                          2.46348628146576e-07, 0.482570902122719, 1958893.73069661, 1958894.21326752
                          ), .Names = c("R0", "r", "infper", "i0", "I0", "S0", "N"))

    expect_equal(summarize.pars(trans.pars(ss0, "unconstrained")),
                 oldval,
                 tolerance=1e-6)
    f1 <- fitsir(bombay2,start=ss0, method="Nelder-Mead", optimizer="optim")
    oldval <- structure(c(1.03206001777913, 0.37953835593663, 0.0844710877771807, 
                          2.59211314281603e-06, 4.59204375483407, 1771539.897656, 1771544.48969976
                          ), .Names = c("R0", "r", "infper", "i0", "I0", "S0", "N"))

    expect_equal(summarize.pars(coef(f1)),
                 oldval,
                 tolerance=1e-6)
})

test_that("philadelphia fits", {
    phila1918a <- with(phila1918, data.frame(times=seq_along(date), count=pim))
    ss0 <- startfun(data=phila1918a, type="death")
    oldval <- structure(c(1.15431358642151, 0.200572186303702, 0.769366826304856, 
                          5.68553945299452e-06, 0.254251231864238, 44718.6741752214, 44718.9284264533
                          ), .Names = c("R0", "r", "infper", "i0", "I0", "S0", "N"))
    
    expect_equal(summarize.pars(trans.pars(ss0, "unconstrained")),
                 oldval,
                 tolerance=1e-6)
    
    suppressWarnings(f1 <- fitsir(phila1918a,start=ss0,type="death", method="BFGS", optimizer="optim"))
    oldval <- structure(c(2.30011752489424, 0.316692031313227, 4.10530545875451, 
                          3.52493040787634e-06, 0.0530463797851002, 15048.8624347797, 15048.9154811595
                          ), .Names = c("R0", "r", "infper", "i0", "I0", "S0", "N"))
    
    ## why are they different??? what happened??
    expect_equal(summarize.pars(coef(f1)),
                 oldval,
                 tolerance=1e-6)
        
})

test_that("harbin fits", {
    ss0 <- startfun(data=harbin, tcol="week", icol="Deaths", type="death")
    oldval <- structure(c(1.93949639324649, 0.781494541260066, 1.20217908589824, 
                          0.000614744499507107, 0.908076582157615, 1476.25289498521, 1477.16097156737
                          ), .Names = c("R0", "r", "infper", "i0", "I0", "S0", "N"))
    
    expect_equal(summarize.pars(trans.pars(ss0, "unconstrained")),
                 oldval,
                 tolerance=1e-6)
    
    suppressWarnings(f1 <- fitsir(harbin,start=c(ss0, phi=5), 
                                  type="death",method="BFGS", 
                                  family="nbinom1",tcol="week",icol="Deaths"))
    oldval <- structure(c(1.86311570920975, 0.797999294064676, 1.08159958991117, 
                          0.000351213159194791, 0.700633492159998, 1994.19469949111, 1994.89533298327
                          ), .Names = c("R0", "r", "infper", "i0", "I0", "S0", "N"))
    expect_equal(summarize.pars(coef(f1)),
                 oldval,
                 tolerance=1e-6)
    
})
