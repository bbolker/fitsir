library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)
library(fitsir)
load("stochsim_data_fix.rda")

file.save <- FALSE

nsim <- length(simlist)

targetpars <- names(fixedpars)

testdist <- c("gaussian", "poisson", "quasipoisson", "nbinom")

level <- 0.95
ll <- (1-level)/2
z.default <- qnorm(1-ll)

confint_summary <- function(x, z=z.default) {
    ss <- summary(x)
    t(data.frame(lwr=ss@summary[1,]-z*ss@summary[2,], 
               upr=ss@summary[1,]+z*ss@summary[2,]))
}

contain <- function(x, interval) interval[1] <= x & interval[2] >= x

filename <- "compare_dist.rda"

if(file.exists(filename)) {
    load(filename)
} else {
    fitList <- vector("list", nsim)
    for(i in 1:nsim) {
        df <- simlist[[i]]
        fitList[[i]] <- Map(fitsir, dist=testdist, MoreArgs=list(data=df, type="incidence", method="BFGS"))
    }
    
    if (file.save) save("fitList", file="compare_dist.rda")
}

if(!exists("resList")) {
    resList <- vector("list", nsim)
    
    for (i in 1:nsim) {
        ## print(i)
        res <- fitList[[i]]
        
        summ <- lapply(res, function(x) {
            conf <- confint_summary(x)
            conf <- as.data.frame(conf)
            coef <- as.data.frame(t(summary(x)@summary[1,]))
            
            ct <- lapply(targetpars, function(name) {
                data.frame(coverage=contain(fixedpars[name], conf[[name]]),
                           width=diff(conf[[name]]),
                           estimate=mean(coef[[name]]))
            })
            
            names(ct) <- targetpars
            
            ct %>% bind_rows(.id="param")
        })
        tmp <- summ %>% bind_rows(.id="dist")
        resList[[i]] <- tmp
        
        if (file.save) save("fitList", "resList", file="compare_dist.rda")
        ## print(tmp)
    }
}

combList <- resList %>% 
    bind_rows(.id="sim") %>%
    group_by(dist, param)

covList <- combList %>%
    mutate(coverage=ifelse(is.na(coverage), FALSE, coverage)) %>%
    summarise(value=sum(coverage)/nsim) %>%
    mutate(key="coverage", 
           lwr=value-z.default*sqrt(value*(1-value)/nsim),
           upr=value+z.default*sqrt(value*(1-value)/nsim)) 

covList %>% 
    select(dist, param, value) %>%
    spread(dist, value) %>% 
    print

mL <- combList %>%
    select(-coverage) %>%
    gather(key, value, -dist, -param, -sim)

g.base <- ggplot(mL, aes(x=param, y=value, group=interaction(dist, param), col = dist)) +
    facet_wrap(~key, scale="free")

g1 <- g.base +
    geom_boxplot() +
    scale_y_log10() +
    theme(axis.title=element_blank())

g2 <- g.base %+% covList +
    geom_point(position=position_dodge((width=0.5))) +
    geom_errorbar(aes(ymax=upr,ymin=lwr), position=position_dodge((width=0.5))) +
    geom_hline(yintercept=level, lty=2) +
    theme(legend.position="none",
          axis.title.x=element_blank())

grid.arrange(g2, g1, nrow = 1, widths=c(1,2.4))

## sample simulation

j <- 10
df <- simlist[[j]]
ff <- fitList[[j]]$gaussian
ff2 <- fitList[[j]]$nbinom

plot(ff2, level=0.95, col.traj="red")
plot(ff, level=0.95, col.traj="blue", col.conf="blue", add=TRUE)
legend(0, 450, c("NB2", "Gaussian"),col=c("red", "blue"), lty=1)

## confidence interval based on mvrnorm
par(mfrow=c(1,2))
predict(ff, level=0.95, method="mvrnorm", debug=TRUE)
predict(ff2, level=0.95, method="mvrnorm", debug=TRUE)
