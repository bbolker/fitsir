g <- SIR.logLik(incidence = FALSE)
tmpf <- function(log.beta,log.gamma,basepars=fpars,data=bombay2,
                 useSSQ=FALSE) {
    pars <- basepars
    pars[c("log.beta","log.gamma")] <- c(log.beta,log.gamma)
    if (useSSQ) return(findSSQ(data,pars)$SSQ)
    return(g(pars,data$count))
}
summarize.pars(fpars)
tmpf2 <- function(log.R0plus1,log.r,basepars=fpars,data=bombay2) {
    pars <- basepars
    pars[c("log.beta","log.gamma")] <- c(log.beta,log.gamma)
    findSSQ(data,pars)$SSQ
}

tmpf(2.43,2.4)
cc <- curve3d(tmpf(x,y),xlim=c(2.4,2.44),
              ylim=c(2.37,2.42),n=c(41,41))
ccssq <- curve3d(tmpf(x,y,useSSQ=TRUE),xlim=c(2.4,2.44),
              ylim=c(2.37,2.42),n=c(41,41))
cc2 <- curve3d(tmpf(x,y),xlim=c(2.42,2.435),
              ylim=c(2.38,2.40),n=c(51,51))
cc3 <- curve3d(tmpf(x,y,useSSQ=TRUE),xlim=c(2.42,2.435),
              ylim=c(2.38,2.40),n=c(51,51))

cc$z[cc$z>1e3] <- NA
cc$z2 <- cc$z

which(cc$z>800,arr.ind=TRUE)
cc$x[9]
cc$y[3]
points(cc$x[9],cc$y[3])
tmpf(cc$x[9],cc$y[3])
plot(count~tvec,data=bombay2)
pp <- pars
pp["log.beta"] <- cc$x[9]
pp["log.gamma"] <- cc$y[3]
ss <- SIR.detsim(bombay2$tvec,trans.pars(pp))
lines(bombay2$tvec,ss)

with(ccssq,image(x,y,log10(z-min(z,na.rm=TRUE)+100),
              xlab="beta",ylab="gamma"))

ss2 <- SIR.detsim(bombay2$tvec,trans.pars(pp),method="rk4",hini=0.001)

with(ccssq,persp(x,y,log10(z-min(z)+100)))
abline(a=0,b=1)
library(rgl)
with(ccssq,persp3d(x,y,log10(z-min(z,na.rm=TRUE)+100),col="gray"))
open3d()
with(cc,persp3d(x,y,log10(z-min(z,na.rm=TRUE)+10),col="blue"))
open3d()
with(cc2,persp3d(x,y,log10(z-min(z,na.rm=TRUE)+0.0001),col="green"))
open3d()
with(cc3,persp3d(x,y,log10(z-min(z,na.rm=TRUE)+10),col="yellow"))
with(cc,persp(x,y,z,ticktype="detailed",zlim=c(96000,1e6)))


