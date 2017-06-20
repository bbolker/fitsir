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
(iniI <- unname(exp(predict(m, data.frame(times=times[1])))))

## ------------------------------------------------------------------------
Qp.alt <- predict(ss,ss.tmax,deriv=2)$y ## second derivative
Ip <- exp(max(predict(ss)$y)) ## I_peak
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
t.diff <- diff(times)
t.diff <- c(t.diff[1], t.diff)
count.orig <- count
count <- count/t.diff

## ------------------------------------------------------------------------
death.pars <- list()
(death.pars <- within(prev.pars,{
    gamma <-  c*(d0 + Ip)/r
    beta <- gamma + r
    N <- beta*gamma/c
    i0 <- iniI/(gamma*N)
}))

plot(times, count, log="y")
lines(times, SIR.detsim(times, unlist(death.pars)[c(3, 4, 2, 1)], type="death"), col="red")

## ----final size----------------------------------------------------------
ss.t3 <- floor(ss.tmax+0.25*ss.tmax)
ss.t4 <- ceiling(ss.tmax+0.75*ss.tmax)
m2 <- lm(log(count)~times,data=subset(ss.data,times<=ss.t4 & times>=ss.t3))
        
times.predict1 <- seq(min(times), ss.t2, by = t.diff[length(t.diff)])
times.predict2 <- seq(ceiling(ss.t3), 3*ss.tmax, by = t.diff[length(t.diff)])
count.predict1 <- exp(predict(m, data.frame(times = times.predict1)))
count.predict2 <- exp(predict(m2, data.frame(times = times.predict2)))
finalsize <- sum(count.predict1) + sum(count.orig[times > ss.t2 & times <= ss.t3]) + sum(count.predict2)

## ------------------------------------------------------------------------
sizefun <- function(beta) {
    R0 <- beta/(beta-r)
    N <- (beta + r)/c
    (finalsize/N - (1 - exp(-R0 * finalsize/N)))^2
}
betavec <- seq(1.1 * r, 50*r, r/10)
sizevec <- sizefun(betavec)

if (all(diff(sizevec) < 0)) {
    beta <- betavec[head(which(sizevec < 1e-4), 1)]
} else {
    beta <- betavec[which.min(sizevec)]
}

## ------------------------------------------------------------------------
inc.pars <- list()
(inc.pars <- within(prev.pars,{
    beta <- beta
    gamma <- beta - r
    N <- (beta + r)/c
    i0 <- iniI/beta/N
}))

plot(times, count, log="y")
lines(times, SIR.detsim(times, unlist(inc.pars)[c(3, 4, 2, 1)], type="incidence"), col="red")

