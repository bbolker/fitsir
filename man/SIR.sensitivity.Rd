% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitsir.R
\name{SIR.sensitivity}
\alias{SIR.sensitivity}
\title{Gradient of negative log likelihood with respect to each parameters}
\usage{
SIR.sensitivity(params, count, times = NULL, model,
  type = c("prevalence", "incidence", "death"), debug = FALSE)
}
\arguments{
\item{params}{parameter vector (log.N, logit.i, log.beta, log.gamma)}

\item{count}{data (epidemic counts for each time period)}

\item{times}{time vector}

\item{model}{log-likelihood model}

\item{type}{type of reported data}

\item{debug}{print debugging output?}
}
\description{
Gradient of negative log likelihood with respect to each parameters
}
