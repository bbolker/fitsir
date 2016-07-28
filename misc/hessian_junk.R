library(fitsir)
library(emdbook)
findSSQ <- fitsir:::findSSQ
findSens <- fitsir:::findSens
fitsir.optim <- fitsir:::fitsir.optim
fpars <- structure(c(2.51906962912771, 2.48802919650103, 14.4215000824118, 
-12.8987594479826), .Names = c("log.beta", "log.gamma", "log.N", 
"logit.i"))
fpars2 <- structure(c(1.40166238626354, 1.30139297147845, 12.1152292130852, 
-10.6845856323438), .Names = c("log.beta", "log.gamma", "log.N", 
"logit.i"))

bombay2 <- setNames(bombay,c("tvec","count"))

g <- fitsir:::SIR.logLik(incidence = FALSE)
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
cc4 <- curve3d(tmpf(x,y),xlim=c(2.38,2.46),
              ylim=c(2.3,2.4),n=c(41,41))

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
with(cc,image(x,y,log10(z-min(z,na.rm=TRUE)+10)))
with(cc3,image(x,y,log10(z-min(z,na.rm=TRUE)+10)))
with(cc,persp(x,y,z,ticktype="detailed",zlim=c(96000,1e6)))
with(cc4,image(x,y,log10(z-min(z,na.rm=TRUE)+10)))

## identifying the appropriate transect via locator(2):
## is there some more automated way to do this??
## strans <- cbind(log.beta=seq(2.42,2.433,length.out=201),
##                 log.gamma=seq(2.386,2.400,length.out=201))
plot(cc4$z[,1])
plot(cc4$z[1,])
plot(cc4$z[,nrow(cc4$z)])
w1 <- which.min(cc4$z[1,])
w2 <- which.min(cc4$z[,nrow(cc4$z)])
## p0 <- c(2.404,2.369)
## p1 <- c(2.4405,2.4066)
p0 <- c(cc4$x[1],cc4$x[w2])
p1 <- c(cc4$y[w1],cc4$y[length(cc4$y)])
ss <- seq(-0.5,1.5,length=101)
strans <- cbind(log.beta=ss*(p1[1]-p0[1])+p0[1],
                log.gamma=ss*(p1[2]-p0[2])+p0[2])

strans.nll <- apply(strans,1,function(p) tmpf(p[1],p[2]))
plot(strans.nll)
## extend the transect?

