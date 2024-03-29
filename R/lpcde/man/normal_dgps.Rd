% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lpbwcde_fns.R
\name{normal_dgps}
\alias{normal_dgps}
\title{Internal function.}
\usage{
normal_dgps(x, v, mean, sd)
}
\arguments{
\item{x}{Scalar, point of evaluation.}

\item{v}{Nonnegative integer, the derivative order (0 indicates cdf, 1 indicates pdf, etc.).}

\item{mean}{Scalar, the mean of the normal distribution.}

\item{sd}{Strictly positive scalar, the standard deviation of the normal distribution.}
}
\description{
Calculates density and higher order derivatives for Gaussian models.
}
\keyword{internal}
