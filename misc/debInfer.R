library(deSolve)
library(deBInfer)
library(fitsir)

f <- fitsir.deBInfer(bombay2,
    R0.args = list(value = 1.007, hypers = list(min = 1, max = 1.01)),
    r.args = list(value = 0.4, hypers = list(min = 0.3, max = 0.5)),
    S.args = list(value = 3e7, hypers = list(min = 2e7, max = 4e7)),
    I.args = list(value = 5, hypers = list(min =4, max = 6)),
    iter = 5000,
    plot = TRUE)

plot(f)

burnin <- 2000
pairs(f, burnin = burnin, scatter=TRUE, trend=TRUE)

post_traj <- post_sim(f, n=500, times=bombay2$tvec, burnin=burnin, output = 'all', prob = 0.95)

## this doesn't work with bombay data for some reason
plot(post_traj, plot.type = "medianHDI", lty = c(1,2), lwd=3, col=c("red","grey20"),
		panel.first=lines(data, col="darkblue", lwd=2))

## not getting a good fit...

plot(bombay2, log = "y")
lines(bombay2$tvec, post_traj[["median"]]$I)
matlines(bombay2$tvec, post_traj[["HDI"]]$I, col = 2, lty = c(2,2))

## TODO: Write a wrapper function for beta/gamma/N/i0 distribution...?

## best fit
f$samples[which.max(f$lpost),]

summary(f)
