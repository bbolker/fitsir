library(ggplot2); theme_set(theme_bw())
library(gridExtra)
library(fitsir)
library(dplyr)
library(tidyr) #gather
library(magrittr) ## %<>%
library(tibble) ## add_column
load("../LHSsim_full.rda")

zero_margin <- theme(panel.margin=grid::unit(0,"lines"))
no_legend <- theme(legend.position = "none")
scale_colour_discrete <- function(...,palette="Dark2") scale_colour_brewer(...,palette=palette)

getsteps <- . %>%
    lapply(function(x){
        df <- x$fitted2
        n <- (df$ll %>% sort %>% diff) > 1
        sum(n)
    }) %>% 
    unlist

steps <- tibble(NM = getsteps(resList), BFGS = getsteps(resList.grad))

steps2 <- steps %>% group_by(NM, BFGS) %>% summarize(rep = length(BFGS))

## For some reason, Nelder-Mead gives us way more local minimas and get stuck in weird places...
## For now, I'm just going to look at BFGS fits...
ggplot(steps2, aes(NM, BFGS, size = rep)) + 
    geom_point() + 
    geom_abline(intercept = 0, slope = 1, col = 2, lty = 2)

i <- steps$BFGS %>% which.max

start <- resList.grad[[i]]$start
fits <- resList.grad[[i]]$fitted2
data <- resList.grad[[i]]$data

fits2 <- fits %>%
    as_data_frame %>%
    mutate(nll = -ll, r = rank(-ll, ties.method = "first"))
    
tmp <- fits2$nll %>%
    round %>%
    unique %>%
    sort

cuts <- (tmp[-1] + tmp[-length(tmp)])/2
cuts <- c(0, cuts, Inf)
cc <- cut(fits2$nll, breaks = cuts)
cc <- (1:length(tmp))[cc]

fits2 %<>% add_column(cc)

example <- fits2 %>% 
    group_by(cc) %>%
    summarize(r = round((min(r) + max(r))/2))

ex.fits <- fits2 %>% 
    filter(r %in% example$r) %>%
    arrange(r)

ex.sim <- ex.fits %>%
    apply(1, function(x){
        tibble(tvec = data$tvec, count = SIR.detsim(data$tvec, trans.pars(x[1:4]), incidence = TRUE))
    }) %>%
    bind_rows(.id = "run")

label <- round(ex.fits$nll, 2)

g.nll <- ggplot(fits2, aes(r, nll)) + 
    geom_line(lty = 2, lwd = 0.3) +
    geom_line(aes(col = factor(cc))) +
    geom_text(data = ex.fits, label = label, vjust = -0.4) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(limits = c(min(tmp)-1, max(tmp)+20)) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    no_legend

g.traj <- ggplot(ex.sim, aes(tvec, count)) +
    geom_line(aes(col = run)) +
    annotate("text", x = 0, y = max(data$count), label = label, hjust = 0.1, vjust = 1.2) +
    geom_point(data = data, size = 0.8) +
    scale_x_continuous(name = "time") +
    facet_wrap(~run, nrow = 1) +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank()) +
    zero_margin +
    no_legend

g.comb <- arrangeGrob(g.nll, g.traj, nrow = 2, heights = c(1, 2.5))

ggsave("nll_steps.pdf", g.comb, width = 8, height = 6, dpi = 600)

summarize.pars2 <- . %>%
    apply(1, function(x){
        x[1:4] %>% summarize.pars %>% t %>% as_tibble
    }) %>%
    bind_rows(.id = "run")

combList <- list(begin = start, end = fits2)

combSum <- Map(summarize.pars2, combList) %>%
    bind_rows(.id = "fit")

ggplot(combSum) +
    geom_point(aes(R0, r, col = fit)) +
    scale_y_log10()

