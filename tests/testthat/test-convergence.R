stopifnot(require("testthat"), require("fitsir"))

context("convergence tests")
test_that("convergence", {
    harbin2 <- setNames(harbin,c("times","count"))
    
    for (dist in c("gaussian", "nbinom", "nbinom1", "poisson")) {
        ss <- startfun(harbin2)
        
        if(dist=="nbinom") ss <- c(ss, ll.k=5)
        
        if(dist=="nbinom1") ss <- c(ss, ll.phi=5)    
        
        f1 <- fitsir(harbin2, type="death", start=ss, dist=dist)
        f2 <- fitsir(harbin2, type="death", method="BFGS", start=coef(f1), dist=dist)
        f3 <- fitsir(harbin2, type="death", start=coef(f2), dist=dist)
        f4 <- fitsir(harbin2, type="death", method="BFGS", start=coef(f3), dist=dist)
        
        expect_equal(coef(f3),
                     coef(f4),
                     tolerance=1e-6)
        
        expect_equal(
            all(eigen(f4@details$hessian)$values>0),
            TRUE
        )
        
        expect_equal(
            f4@details$convergence,
            0
        )
    }
})
