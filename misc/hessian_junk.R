library(devtools)
load_all("..")

bombay2 <- setNames(bombay, c("tvec", "count"))
fpars <- coef(f1 <- fitsir(bombay2, start = startfun()))

##more stable than fpars?
fpars2 <- coef(f2 <- fitsir(bombay2, start = fpars))

g <- SIR.logLik(incidence = FALSE, method = "rk4", hini = 0.1)
tmpf <- function(log.beta,log.gamma,basepars=fpars,data=bombay2,
                 useSSQ=FALSE) {
    pars <- basepars
    pars[c("log.beta","log.gamma")] <- c(log.beta,log.gamma)
    if (useSSQ) return(findSSQ(data,pars)$SSQ)
    return(g(pars,data$count))
}

library(emdbook)

fpars2

cc <- curve3d(tmpf(x,y, basepars = fpars2),xlim=c(2.518,2.53),
              ylim=c(2.487,2.49),n=c(30,30))

yminvals <- apply(cc$z,1,function(x) cc$y[which.min(x)])
xminvals <- apply(cc$z,1,function(x) cc$x[which.min(x)])

plot(xminvals,yminvals)

library(rgl)
open3d()
with(cc,persp3d(x,y,log10(z-min(z,na.rm=TRUE)+0.001),col="blue"))

min(cc$z)

which(cc$z == min(cc$z), arr.ind = TRUE)

m.x <- cc$x[19]
d.x <- cc$x[19] - cc$x[17]

m.y <- cc$y[19]
d.y <- cc$y[19] - cc$y[17]

cc.zoom <- curve3d(tmpf(x,y, basepars = fpars2),xlim=c(m.x - d.x, m.x + d.x),
                   ylim=c(m.y-d.y,m.y+d.y),n=c(21,21))

open3d()
with(cc.zoom,persp3d(x,y,log10(z-min(z,na.rm=TRUE)+0.01),col="red"))

with(cc.zoom,persp(x,y,log10(z-min(z)+100)))
with(cc.zoom,image(x,y,log10(z-min(z,na.rm=TRUE)),
               xlab="beta",ylab="gamma"))

cc2 <- cc
cc2$z[cc2$z > 167.9] <- NA

with(cc,image(x,y,log10(z-min(z,na.rm=TRUE)+1e-5),
                   xlab="beta",ylab="gamma"))


##very different....
fpars3 <- coef(f1 <- fitsir(bombay2, start = fpars2))