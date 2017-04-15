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
z <- qnorm(1-ll)

confint_summary <- function(x, z=z) {
    ss <- summary(x)
    t(data.frame(lwr=ss@summary[1,]-z*ss@summary[2,], 
               upr=ss@summary[1,]+z*ss@summary[2,]))
}

contain <- function(x, interval) interval[1] <= x & interval[2] >= x

fitList <- vector("list", nsim)
resList <- vector("list", nsim)

filename <- "compare_dist.rda"

if(file.exists(filename)) {
    load(filename)
} else {
    for(i in 1:nsim) {
        df <- simlist[[i]]
        fitList[[i]] <- Map(fitsir, dist=testdist, MoreArgs=list(data=df, type="incidence", method="BFGS"))
    }
}

if(file.exists(filename)) {
    load(filename)
} else {
    for (i in 1:nsim) {
        ## print(i)
        res <- fitList[[i]]
        
        conf <- lapply(res, confint_summary)
        
        summ <- lapply(conf, function(x) {
            x <- as.data.frame(x)
            ct <- lapply(targetpars, function(name) {
                data.frame(coverage=contain(fixedpars[name], x[[name]]),
                           width=diff(x[[name]]),
                           estimate=mean(x[[name]]))
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
           lwr=value-z*sqrt(value*(1-value)/nsim),
           upr=value+z*sqrt(value*(1-value)/nsim)) 

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
