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

