library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(optimx)
library(fitsir)
source("optimx.R")
load("stochsim_data_fix.rda")

file.save <- TRUE

nsim <- length(simlist)

targetpars <- names(fixedpars)

filename <- "compare_method.rda"

if(file.exists(filename)) {
    load(filename)
} else {
    fitList <- vector("list", nsim)
    for(i in 1:nsim) {
        cat(i)
        df <- simlist[[i]]
        ss <- c(startfun(df, type="incidence"), log.dsp=3)
        suppressWarnings(fitList[[i]] <- try(fitsirx(df, start=ss, type="incidence", dist="nbinom")))
        if (file.save) save("fitList", file=filename)
    }
    
    if (file.save) save("fitList", file=filename)
}

(is_err <- fitList %>%
    lapply(inherits, "try-error") %>%
    unlist %>%
    which)
    
is_null <- fitList %>%
    lapply(is.null) %>%
    unlist %>%
    which

fitList2 <- fitList[-c(is_err, is_null)]

goodfits <- lapply(fitList2, function(x){
    ss <- summary(x, order="value")
    ss[which(min(ss$value)*1.01 > ss$value),]
})

goodmethods <- goodfits %>%
    lapply(rownames) %>%
    unlist

## So BFGS performs fairly well if we start on the right parameters...
## Notice that nlminb performs well on all simulations except 1
goodmethods %>%
    unlist %>%
    table

lapply(goodfits, function(x){
    summary <- apply(x, 1, summarize.pars)
    data.frame(method=colnames(summary), t(summary))
}) %>% bind_rows %>%
    gather(key, value,-method) %>%
    filter(key %in% c("R0", "infper", "N", "I0")) %>%
    ggplot(aes(method, value, col=method)) +
        geom_violin() +
        scale_y_log10() +
        facet_wrap(~key, scale="free")



