library(dplyr)
library(ggplot2)
library(fitsir)
load("stochsim_data_fix.rda")

save <- FALSE

nsim <- length(simlist)

targetpars <- names(fixedpars)

testdist <- c("gaussian", "poisson", "quasipoisson", "nbinom")

confint_summary <- function(x, level=0.95) {
    ss <- summary(x)
    ll <- (1-level)/2
    z <- qnorm(1-ll)
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
        
        if (save) save("fitList", "resList", file="compare_dist.rda")
        ## print(tmp)
    }
}

combList <- resList %>% 
    bind_rows(.id="sim") %>%
    group_by(dist, param)

covList <- combList %>%
    mutate(coverage=ifelse(is.na(coverage), FALSE, coverage)) %>%
    summarise(coverage=sum(coverage)/nsim)

covList %>% 
    spread(dist, coverage) %>% 
    print

mL <- combList %>%
    select(-coverage) %>%
    gather(key, value, -dist, -param, -sim)

ggplot(mL) +
    geom_boxplot(aes(x=param, y=value, group=interaction(dist, param), col=dist)) +
    scale_y_log10() +
    facet_wrap(~key)
