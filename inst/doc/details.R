## ----opts,echo=FALSE-----------------------------------------------------
library("knitr")
opts_chunk$set(fig.width=5,fig.height=5,tidy=FALSE,message=FALSE,error=FALSE,warning=FALSE)

## ----spline--------------------------------------------------------------
library(fitsir)
times <- seq_along(phila1918$date)
count <- phila1918$pim
ss <- fitsir:::smooth.spline2(times, count, itmax = 100, relpeakcrit = 0.1)

## ----little_r------------------------------------------------------------
ss.data <- data.frame(times = times, count = exp(predict(ss)$y))
ss.tmax <- ss.data$times[which.max(ss.data$count)]
ss.t1 <- min(times)+0.25*(ss.tmax-min(times))
ss.t2 <- min(times)+0.75*(ss.tmax-min(times))
        
m <- lm(log(count)~times,data=subset(ss.data,times<=ss.t2 & times>=ss.t1))
(r <- unname(coef(m)[2])) ##beta - gamma

## ----little_r_plot, echo=FALSE-------------------------------------------
plot(times, count, log="y")
lines(ss.data)
lines(11:31, exp(predict(m)), col="red")

## ----ini_size------------------------------------------------------------
(iniI <- unname(exp(coef(m)[1])))

## ------------------------------------------------------------------------
Qp.alt <- predict(ss,ss.tmax,deriv=2)$y
Ip <- exp(max(predict(ss,times)$y))
(c <- -Qp.alt/Ip)

## ----integral------------------------------------------------------------
t.diff <- diff(times)
t.diff <- c(t.diff[1], t.diff)
ss.int <- transform(ss.data, int = count * t.diff)
ss.int <- ss.int[times<ss.tmax,]

d0 <- sum(ss.int[,3])

## ----loop----------------------------------------------------------------
while(r - c * d0 < 0){
    ss.int <- ss.int[-nrow(ss.int),]
    d0 <- sum(ss.int[,3])
}

## ------------------------------------------------------------------------
prev.pars <- list()
(prev.pars <- within(prev.pars,{
    gamma <- c * Ip/(r - c * d0)
    beta <- gamma + r
    N <- beta*gamma/c
    i0 <- iniI/N
}))

plot(times, count, log="y")
lines(times, SIR.detsim(times, unlist(prev.pars)[c(3, 4, 2, 1)]), col="red")

## ------------------------------------------------------------------------
death.pars <- list()
(death.pars <- within(prev.pars,{
    gamma <-  c*(d0 + Ip)/r
    beta <- gamma + r
    N <- beta*gamma/c
    i0 <- iniI/N
}))

plot(times, count, log="y")
lines(times, SIR.detsim(times, unlist(death.pars)[c(3, 4, 2, 1)], type="death"), col="red")

