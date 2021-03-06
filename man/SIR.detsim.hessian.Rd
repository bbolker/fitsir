% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hessian.R
\name{SIR.detsim.hessian}
\alias{SIR.detsim.hessian}
\title{integrate second order sensitivities}
\usage{
SIR.detsim.hessian(t, params, type = c("prevalence", "incidence",
  "death"))
}
\arguments{
\item{t}{time vector}

\item{params}{parameter vector (beta, gamma, N, i0)}

\item{type}{type of count data}
}
\description{
integrate second order sensitivities
}
