stopifnot(require("testthat"), require("fitsir"))

context("convergence tests")
test_that("convergence", {
    harbin2 <- setNames(harbin,c("times","count"))
    
    for (dist in c("gaussian", "nbinom", "nbinom1", "poisson")) {
        ss <- startfun(harbin2)
        
        if(dist=="nbinom") ss <- c(ss, ll.k=5)
        
        if(dist=="nbinom1") ss <- c(ss, ll.phi=5)    
        
        f1 <- fitsir(harbin2, type="death", start=ss, dist=dist)
        
        expect_equal(
            all(eigen(f1@details$hessian)$values>0),
            TRUE
        )
        
        expect_equal(
            f1@details$convergence,
            0
        )
    }
})
