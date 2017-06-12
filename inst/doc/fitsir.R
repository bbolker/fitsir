## ----opts,echo=FALSE-----------------------------------------------------
library("knitr")
opts_chunk$set(fig.width=5,fig.height=5,tidy=FALSE,message=FALSE,error=FALSE,warning=FALSE)

## ----pkgs,message=FALSE--------------------------------------------------
library(fitsir)
library(dplyr)
library(ggplot2); theme_set(theme_bw())

## ----harbin_head---------------------------------------------------------
head(harbin)

## ----harbin_fit----------------------------------------------------------
dietz_harbin <- c(x0=2985,rzero=2,gamma=7/11)
dietz_lpars <- with(as.list(dietz_harbin),
      c(log.beta=log(rzero*gamma),
        log.gamma=log(gamma),
        log.N=log(x0),
        logit.i=qlogis(1e-3)))
(ff <- fitsir(harbin, start=dietz_lpars, type="death", 
              tcol="week", icol="Deaths", method="BFGS"))

## ----harbin_plot, echo = TRUE, message = FALSE, fig.height=4-------------
plot(ff, main="SIR vs. KM comparison")
times <- harbin$week
dpKM <- with(as.list(dietz_harbin), 
        {
           c1 <- sqrt((rzero-1)^2+2*rzero^2/x0)
           c2 <- atanh((rzero-1)/c1)
           gamma*x0/(2*rzero^2)*c1*
               (1/cosh(c1*gamma*times/2-c2))^2
        })
lines(times,dpKM, col = 2)
legend(x=2, y=275, legend=c("SIR","Dietz"), col=c("black", "red"), lty = 1)

## ----harbin_summary------------------------------------------------------
summary(ff)

## ----harbin_test---------------------------------------------------------
dispersion(ff)

## ----harbin2-------------------------------------------------------------
harbin2 <- setNames(harbin, c("times", "count"))

## ----harbin_dispersion---------------------------------------------------
ff2 <- fitsir(harbin2, dist="quasipoisson", type="death", method="BFGS", start=dietz_lpars)
ff3 <- fitsir(harbin2, dist="nbinom", type="death", method="BFGS", start=c(dietz_lpars, ll.k=5))
ff4 <- fitsir(harbin2, dist="nbinom1", type="death", method="BFGS", start=c(dietz_lpars, ll.phi=5))

## ----harbin_dispersion_plot, fig.width=8---------------------------------
plot(ff2, level=0.95, col.traj="green", col.conf="green", log="y", main="Comparison of three error functions")
plot(ff3, level=0.95, add=TRUE, col.traj="blue", col.conf="blue")
plot(ff4, level=0.95, add=TRUE, col.traj="red", col.conf="red")
legend(x=2, y=275, legend=c("Qausipoisson","NB2", "NB1"), col=c("green", "blue", "red"), lty = 1)

## ----harbin_logLik-------------------------------------------------------
hfits <- list(QP=ff2, NB2=ff3, NB1=ff4)
lapply(hfits, AIC)

## ----mean-variance, echo=FALSE-------------------------------------------
mvrel <- function(fit, data) {
    mean <- SIR.detsim(data$times, coef(fit, "trans"), type="death")
    data.frame(
        mean=mean,
        var=(data$count-mean)^2
    )
}
level <- seq(0, 300, by = 25)

mvfun <- . %>%
    mvrel(harbin2) %>%
    mutate(group=cut(mean, breaks=level)) %>%
    group_by(group) %>%
    summarise(mean2 = mean(mean), var2=mean(var), n=length(var))

mvtot <- hfits %>%
    lapply(mvfun) %>%
    bind_rows(.id="dist") %>%
    filter(dist != "QP")

## `dispersion` returns raw dispersion parameters that have not been transformed
nb1k <- dispersion(ff4, "nbinom1")
nb2k <- dispersion(ff3, "nbinom")

ggplot(mvtot, aes(mean2, var2)) + 
    geom_point(aes(size=n, col=dist), pch=1) +
    scale_x_continuous(lim=c(0, 250), expand = c(0,0), name="mean") +
    scale_y_continuous(lim=c(0, 600), expand = c(0, 0), name="variance") +
    scale_size_continuous(range = c(5, 20), guide=FALSE) +
    geom_abline(intercept=0,slope=nb1k+1) + 
    stat_function(fun=function(x) x + x^2/nb2k, linetype=2)

## ----harbin_nb1sum-------------------------------------------------------
summary(ff4)

## ----harbin_acf----------------------------------------------------------
acf(residuals(ff4, "raw"))

## ----phila_head----------------------------------------------------------
head(phila1918)

## ----phila_newdata-------------------------------------------------------
phila1918a <- with(phila1918, data.frame(times=seq_along(date), count=pim)) 
plot(phila1918a, log="y")

## ----phila_fit-----------------------------------------------------------
(pfit <- fitsir(phila1918a, 
    start=c(log.beta=log(0.12), log.gamma=log(0.09), log.N=log(10000), logit.i=qlogis(0.01)),  
    type="death")
)
plot(pfit)

## ----phila_start---------------------------------------------------------
(pstart <- startfun(phila1918a, type="death"))
plot(phila1918a, log="y")
lines(SIR.detsim(phila1918a$times, trans.pars(pstart), type="death"))

## ----phila_start2--------------------------------------------------------
(pfit2 <- fitsir(phila1918a, type="death", start=pstart))

## ----phila_param---------------------------------------------------------
ppars <- as.data.frame(t(sapply(plist, coef)))
col <- c("black", "blue", "red")
ccol <- col[cut(pLik, breaks=c(-900, -600, -520, -500))]

plot(ppars, col=ccol)

## ----phila_traj----------------------------------------------------------
ppred <- plist %>%
    lapply(predict) %>%
    bind_rows(.id="sim") 

ppred$logLik <- rep(pLik, each=length(phila1918a$times))

ggplot(ppred) +
    geom_line(aes(times, mean, group=sim, col=logLik)) +
    geom_point(data=phila1918a, aes(times, count))

## ----phila_best----------------------------------------------------------
pbest <- plist[[which.max(pLik)]]
summary(pbest) ## Nelder-Mead being unstable
pbest2 <- fitsir(phila1918a, type="death", start=coef(pbest), method="BFGS")
summary(pbest2)

## ----bombay_set----------------------------------------------------------
bombay2 <- setNames(bombay, c("times", "count"))
bbstart <- startfun(bombay2, type="death")
bb <- fitsir(bombay2, type="death", dist="nbinom", start=c(bbstart, ll.k=5))

## ----bombay_lm_traj------------------------------------------------------
plot(bombay2)
l <- lapply(blist[goodfits], function(x) plot(x, add=TRUE))

## ----bombay_acf----------------------------------------------------------
acf(residuals(bb, "raw"))

