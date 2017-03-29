##' Fit a spline to an epidemic data
##' @param count data (epidemic counts for each time period)
##' @param times time vector
##' @param maxit maximum number of iterations
##' @param relpeakcrit critical relative peak value to test for single peak
smooth.spline2 <- function(times, count, itmax=100,relpeakcrit=0.1){
    single_peak <- FALSE
    it <- 1
    spar <- 0.5
    while (!single_peak && it<itmax) {
        ss <- smooth.spline(times,log(count),spar=spar)
        dd <- predict(ss,deriv=1)$y
        ## change in sign of first derivative
        dds <- diff(sign(dd))
        spar <- if (spar<0.8) spar+0.05 else (1+spar)/2
        it <- it+1
        ncrit <- sum(dds<0)
        peakvals <- count[dds<0]
        relpeakvals <- peakvals[-1]/peakvals[1]
        single_peak <- ncrit==1 ||
            all(relpeakvals<relpeakcrit)
    }
    if (it==itmax) {
        ## try harder?
        stop("couldn't smooth enough")
    }
    return(ss)
}

##' Starting function
##' @param data data frame with columns \code{times} and \code{count}
##' @param type epidemic data type
##' @param tcol column name for time variable
##' @param icol column name for count variable
##' @param log.beta log of per capita transmission rate
##' @param log.gamma log of recovery/removal rate
##' @param log.N log of population size
##' @param logit.i logit of initial proportion infectious
##' @param itmax maximum number of iterations in \code{\link{smooth.spline2}}
##' @param relpeakcrit critical relative peak value used in \code{\link{smooth.spline2}}
##' @export
startfun <- function(data = NULL,
                     type = c("prevalence", "incidence", "death"),
                     tcol = "times", icol = "count",
                     log.beta=log(0.12),log.gamma=log(0.09),
                     log.N=log(10000),logit.i=qlogis(0.01),
                     itmax=100,relpeakcrit=0.1) {
    if (!is.null(data)) { ## auto start
        type <- match.arg(type)
        
        times <- data[[tcol]]
        count <- data[[icol]]
        ## for smooth.spline(log(count)) ...
        if (any(count<=0)) {
            count <- pmax(count,min(count[count>0])/2)
        }
        
        t.diff <- diff(times)
        t.diff <- c(t.diff, t.diff[length(t.diff)])
        
        if (type != "prevalence") {
            ## this should approximately be equal to the rates
            ## i.e., beta * S * I/N for incidence and gamma * I for death
            count.orig <- count
            count <- count/t.diff
        }
        ## smooth data; start with smoothing par 0.5, try
        ## to increase it until there is a single critical point ...
        ## (check that second deriv is negative???)
        ss <- smooth.spline2(times, count, itmax = itmax, relpeakcrit = relpeakcrit)
        
        ss.data <- data.frame(times = times, count = exp(predict(ss)$y))
        ss.tmax <- ss.data$times[which.max(ss.data$count)]
        ss.t1 <- min(times)+0.25*(ss.tmax-min(times))
        ss.t2 <- min(times)+0.75*(ss.tmax-min(times))
        
        m <- lm(log(count)~times,data=subset(ss.data,times<=ss.t2 & times>=ss.t1))
        r <- unname(coef(m)[2]) ##beta - gamma
        
        ## interpolate starting value based on linear regression
        iniI <- unname(exp(coef(m)[1]))
        
        ## curvature of spline at max
        ## using quadratic fit:
        ## t.sub <- (max(times) - ss.tmax)/2
        ## m4 <- lm(log(count)~poly(times,2,raw = TRUE), data = subset(ss.data, times> ss.tmax - t.sub & times < ss.tmax + t.sub))
        ## Qp.alt <- unname(2*coef(m4)[3])
        Qp.alt <- predict(ss,ss.tmax,deriv=2)$y
        if(Qp.alt > 0){
            stop("second derivative larger than 0")
        }
        Ip <- exp(max(predict(ss,times)$y))
        c <- -Qp.alt/Ip
        
        if (type %in% c("prevalence", "death")) {
            ss.int <- transform(ss.data, int = count * t.diff)
            ss.int <- ss.int[times<ss.tmax, ]
            
            d0 <- sum(ss.int[,3]) - iniI
            
            if (type == "prevalence") {
                while(r - c * d0 < 0){
                    ss.int <- ss.int[-nrow(ss.int),]
                    d0 <- sum(ss.int[,3]) - iniI
                }
                
                gamma <- c * Ip/(r - c * d0)
            } else {
                gamma <- c*(d0 + Ip)/r
            }
            
            beta <- gamma + r
            N <- beta*gamma/c
            i0 <- iniI/N
            
        } else if (type == "incidence") {
            ss.t3 <- floor(ss.tmax+0.25*ss.tmax)
            ss.t4 <- ceiling(ss.tmax+0.75*ss.tmax)
            m2 <- lm(log(count)~times,data=subset(ss.data,times<=ss.t4 & times>=ss.t3))
            
            times.predict1 <- seq(min(times), ss.t2, by = t.diff[length(t.diff)])
            times.predict2 <- seq(ceiling(ss.t3), 3*ss.tmax, by = t.diff[length(t.diff)])
            count.predict1 <- exp(predict(m, data.frame(times = times.predict1)))
            count.predict2 <- exp(predict(m2, data.frame(times = times.predict2)))
            finalsize <- sum(count.predict1) + sum(count.orig[times > ss.t2 & times <= ss.t3]) + sum(count.predict2)
            
            sizefun <- function(beta) {
                R0 <- beta/(beta-r)
                N <- (beta + r)/c
                (finalsize/N - (1 - exp(-R0 * finalsize/N)))^2
            }
            betavec <- seq(1.1 * r, 10*r, r/10)
            sizevec <- sizefun(betavec)
            
            ## FIXME: come up with a better method...
            if (all(diff(sizevec) < 0)) {
                beta <- betavec[head(which(sizevec < 1e-4), 1)]
            } else {
                beta <- betavec[which.min(sizevec)]
            }
            
            gamma <- beta - r
            N <- (beta + r)/c
            i0 <- iniI/beta/N
        }
        
        x <- c(
            log.beta = log(beta),
            log.gamma = log(gamma),
            log.N = log(N),
            logit.i = qlogis(i0)
        )
        
        x
    } else {
        c(log.beta=log.beta,log.gamma=log.gamma,log.N=log.N,
          logit.i=logit.i)
    }
}
