library(fitsir)
library(emdbook)
library(rgl)
library(ggplot2); theme_set(theme_bw())
library(reshape2)
library(dplyr)
library(viridis)

findSSQ <- fitsir:::findSSQ
findSens <- fitsir:::findSens
fitsir.optim <- fitsir:::fitsir.optim

bombay2 <- setNames(bombay,c("tvec","count"))

g <- fitsir:::SIR.logLik(incidence = FALSE)
tmpf <- function(log.beta,log.gamma,basepars=fpars,data=bombay2,
                 useSSQ=FALSE) {
    pars <- basepars
    pars[c("log.beta","log.gamma")] <- c(log.beta,log.gamma)
    if (useSSQ) return(findSSQ(data,pars)$SSQ)
    return(g(pars,data$count))
}

fn <- "hessian_explore.rda"
if (file.exists(fn)) {
    load(fn)
} else {
    fit1 <- fitsir(bombay2,start=startfun(auto=TRUE,data=bombay2))
    fpars1 <- coef(fit1)

    ## 5% slice 
    cc1a <- curve3d(tmpf(x,y,basepars=fpars1),
                  xlim=fpars1["log.beta"]*c(0.95,1.05),
                  ylim=fpars1["log.gamma"]*c(0.95,1.05),
                  n=c(61,61))
    ## 1% slice
    cc1b <- curve3d(tmpf(x,y,basepars=fpars1),
                    xlim=fpars1["log.beta"]*c(0.99,1.01),
                    ylim=fpars1["log.gamma"]*c(0.99,1.01),
                    n=c(61,61))

    ## 0.1% slice
    cc1c <- curve3d(tmpf(x,y,basepars=fpars1),
                    xlim=fpars1["log.beta"]*c(0.999,1.001),
                    ylim=fpars1["log.gamma"]*c(0.999,1.001),
                    n=c(61,61))

    save("fit1","fpars1","cc1a","cc1b","cc1c",file=fn)
}


## appears to be really a local minimum ...
all(eigen(findHess(bombay2,fpars1))$values>0)

## confint(fit1,method="quad") ## fails because mle2 uses
## (bad) numDeriv Hessian

with(cc1b,image(x,y,log10(z-min(z)),xlab="log.beta",ylab="log.gamma"))
abline(a=0,b=1)
text(3.65,3.65,expression(R[0]==1))
## with(cc1a,persp3d(x,y,z,col="gray"))
xmins1a <- apply(cc1a$z,1,min)
xmins1b <- apply(cc1b$z,1,min)
plot(cc1a$x,xmins1a,type="b")
lines(cc1b$x,xmins1b,type="b",col=2)

ymins1a <- apply(cc1a$z,1,min)
ymins1b <- apply(cc1b$z,2,min)
plot(cc1a$y,ymins1a,type="b")
lines(cc1b$y,ymins1b,type="b",col=2)

## transform to R0/r space, plot ...
library(reshape2)

cc.cur <- cc1c
dimnames(cc.cur$z) <- list(log.beta=cc.cur$x,log.gamma=cc.cur$y)
cc.cur.df <- melt(cc.cur$z) %>%
    transmute(nll=value,log.R0=log.beta-log.gamma,
              q=log.beta+log.gamma)

## "q" doesn't really have an obvious epidemiological meaning
## (log(beta*gamma)) - it just happens to be orthogonal to log(R0)

(gg1 <- ggplot(cc.cur.df,aes(log.R0,q,colour=nll))+
    geom_point(size=2)+
    scale_color_viridis())

betw <- function(x,lims) {
    lims[1]<x & x<lims[2]
}
R0lims <- c(0.009,0.0105)
qlims <- c(7.325,7.328)
gg1 %+% subset(cc.cur.df,betw(log.R0,R0lims) & betw(q,qlims)) +
    geom_point(size=15,alpha=0.4)+
    annotate(x=fpars1["log.beta"]-fpars1["log.gamma"],
             y=fpars1["log.beta"]+fpars1["log.gamma"],
             colour="red",pch=1,size=18,
             geom="point")+
    stat_summary_2d(geom="contour",aes(z=nll))

ggplot(cc.cur.df,aes(log.R0,q,z=nll))+
    stat_summary_2d()+scale_fill_viridis()
## would like to add a contour at min(nll)+1.92, but can't quite
## get it done ... stat_summary_2d(geom="contour") gives weird results


summary(subset(cc.cur.df,nll<min(nll)+1.92))





