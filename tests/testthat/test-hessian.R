stopifnot(require("testthat"), require("fitsir"))

context("hessian tests")
test_that("bombay fits", {
    harbin2 <- setNames(harbin, c("times", "count"))
    
    for(type in c("prevalence", "incidence", "death")) {
        print(type)
        ss <- startfun(data=harbin2, type=type)
        
        for (dist in c("gaussian", "poisson", "nbinom", "nbinom1")) {
            print(dist)
            if (dist=="nbinom") {
                ss2 <- c(ss, ll.k=3)
            } else if (dist == "nbinom1") {
                ss2 <- c(ss, ll.phi=3)
            } else {
                ss2 <- ss
            }
            suppressWarnings(ff <- fitsir(harbin2, start=ss2, method="BFGS", dist=dist, type=type))
            hess <- SIR.hessian(harbin2, coef(ff), dist=dist, type=type)
            all.equal(
                hess,
                ff@details$hessian
            )
        }
    }
    
})
