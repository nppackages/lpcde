% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lpbwcde.R
\name{lpbwcde}
\alias{lpbwcde}
\title{Data-driven bandwidth selection for local polynomial conditional density estimators}
\usage{
lpbwcde(
  y_data,
  x_data,
  x,
  y_grid = NULL,
  p = NULL,
  q = NULL,
  grid_spacing = "",
  ng = NULL,
  mu = NULL,
  nu = NULL,
  kernel_type = c("epanechnikov", "triangular", "uniform"),
  bw_type = c("mse-rot", "imse-rot"),
  regularize = NULL
)
}
\arguments{
\item{y_data}{Numeric matrix/data frame, the raw data of independent.}

\item{x_data}{Numeric matrix/data frame, the raw data of covariates.}

\item{x}{Numeric, specifies the evaluation point in the x-direction. Default is median of the dataset.}

\item{y_grid}{Numeric, specifies the grid of evaluation points. When set to default, grid points
will be chosen as 0.05-0.95 percentiles of the data, with a step size of 0.05.}

\item{p}{Nonnegative integer, specifies the order of the local polynomial for \code{Y} used to
construct point estimates. (Default is \code{2}.)}

\item{q}{Nonnegative integer, specifies the order of the local polynomial for \code{X} used to
construct point estimates. (Default is \code{1}.)}

\item{grid_spacing}{String, If equal to "quantile" will generate quantile-spaced grid evaluation points, otherwise will generate equally spaced points.}

\item{ng}{Int, number of grid points to be used in generating bandwidth estimates.}

\item{mu}{Nonnegative integer, specifies the derivative with respect to \code{Y} of the
distribution function to be estimated. \code{0} for the distribution function,
\code{1} (default) for the density funtion, etc.}

\item{nu}{Nonnegative integer, specifies the derivative with respect to \code{X} of the
distribution function to be estimated.}

\item{kernel_type}{String, specifies the kernel function, should be one of
\code{"triangular"}, \code{"uniform"} or \code{"epanechnikov"}.}

\item{bw_type}{String, specifies the method for data-driven bandwidth selection. This option will be
ignored if \code{bw} is provided. Implementable with \code{"mse-rot"} (default, mean squared error-optimal
bandwidth selected for each grid point)}

\item{regularize}{Boolean (default TRUE). Option to regularize bandwidth selection to have atleast
20+max(p, q)+1 datapoints when evaluating the estimator.}
}
\value{
\item{BW}{A matrix containing (1) \code{y_grid} (grid point), (2) \code{bw} (bandwidth)}
\item{opt}{A list containing options passed to the function.}
}
\description{
\code{\link{lpbwcde}} implements the bandwidth selection methods for local
polynomial based conditionaldensity (and derivatives) estimation proposed and studied
in \insertCite{bernoulli}{lpcde}.

Companion command: \code{\link{lpcde}} for estimation and robust bias-corrected inference.

Related \code{Stata} and \code{R} packages useful for nonparametric estimation and inference are
available at \url{https://nppackages.github.io/}.
}
\examples{
# Generate a random sample
set.seed(42);
x_data = rnorm(2000)
y_data = rnorm(2000, mean=x_data)
x = 0

# Construct bandwidth
bw1 <- lpbwcde(y_data = y_data, x_data = x_data, x=x, bw_type = "mse-rot")
summary(bw1)

# Display bandwidths for a subset of y_grid points
summary(bw1, y_grid=bw1$BW[2:5, "y_grid"])

}
\references{
\insertRef{bernoulli}{lpcde}
}
\seealso{
Supported methods: \code{\link{coef.lpbwcde}},
\code{\link{print.lpbwcde}}, \code{\link{summary.lpbwcde}}.
}
\author{
Matias D. Cattaneo, Princeton University. \email{cattaneo@princeton.edu}.

Rajita Chandak (maintainer), Princeton University. \email{rchandak@princeton.edu}.

Michael Jansson, University of California Berkeley. \email{mjansson@econ.berkeley.edu}.

Xinwei Ma, University of California San Diego. \email{x1ma@ucsd.edu}.
}
