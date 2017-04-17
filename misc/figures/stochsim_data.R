## generate data by stochastic simulations
source("../stochsim_funs.R")

nsim <- 500

epi_range <- data.frame(min = c(2, 1, 1e4, 10),
                        max = c(4, 10, 1e5, 100),
                        row.names = c("R0", "infper", "N", "I0"))

ltab <- as.data.frame(
    apply(epi_range, 1,
        function(x) exp(seq(log(x[1]), log(x[2]), length = nsim))
    )
)

set.seed(101)
ltab[] <- lapply(ltab, sample)

getpars <- function(x) with(as.list(x), 
    c(log.beta=log(R0/infper),log.gamma=log(1/infper),
      log.N=log(N),logit.i=qlogis(I0/N))
)

truepars <- t(apply(ltab, 1, getpars))

simpars <- as.data.frame(t(apply(truepars, 1, trans.pars)))
simpars <- lapply(split(simpars,  1:nsim), unlist)

lhslist <- Map(simfun, simpars, MoreArgs=list(seed=101))

save("ltab", "lhslist", file="stochsim_data_lhs.rda")

fixedpars <- c(R0=2, infper=11, N=1e4, I0=100)
simlist <- replicate(nsim, simfun(pars=trans.pars(getpars(fixedpars)), dt=2), simplify = FALSE)

save("fixedpars", "simlist", file="stochsim_data_fix.rda")
