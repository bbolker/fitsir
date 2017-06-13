## ----opts,message=FALSE,echo=FALSE---------------------------------------
library("knitr")
opts_chunk$set(fig.width=4,fig.height=4,error=FALSE)
knit_hooks$set(basefig=function(before, options, envir) {
                   if (before) {
                       ## tweak graphical settings for base figures
                       par(bty="l",las=1)
                   } else { }
               })

## ----pkgs,message=FALSE--------------------------------------------------
library(fitsir)
library(deSolve)
if (packageVersion("odeintr")<"1.5")
    stop("need devel version: try devtools::install_github('thk686/odeintr')")
library(odeintr)
library(Rcpp)
library(ggplot2); theme_set(theme_bw())
library(microbenchmark)

## ----sirgrad0------------------------------------------------------------
SIR.grad

## ----sirgrad2------------------------------------------------------------
SIR.grad2 <- function(t, y, params) {
    list(c(-params[1]*exp(y[2])*y[1]/params[3],
            params[1]*y[1]/params[3]-params[2]))
}
## without division by N
SIR.grad4 <- function(t, y, params) {
    list(c(-params[1]*exp(y[2])*y[1],
            params[1]*y[1]-params[2]))
}
## define Jacobian, maybe useful later ...
jacfunc <- function(t, y, params) {
    matrix(c(-exp(y[2])*y[1],0,
             y[1],-1),nrow=2,byrow=TRUE)
}

## ----compCpp-------------------------------------------------------------
s <- sourceCpp("sirgrad.cpp")

