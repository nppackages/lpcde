% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lpbwcde_methods.R
\name{summary.lpbwcde}
\alias{summary.lpbwcde}
\title{Summary method for local polynomial conditional density bandwidth selection}
\usage{
\method{summary}{lpbwcde}(object, ...)
}
\arguments{
\item{object}{Class "lpbwcde" object, obtained by calling \code{\link{lpbwcde}}.}

\item{...}{Additional options, including (i) \code{y_grid} specifies a subset of y_grid points
to display the bandwidth; (ii) \code{gridIndex} specifies the indices of y_grid points
to display the bandwidth.}
}
\value{
\item{Display output}{A list of specified options and a matrix of grid points, bandwidth, and effective sample size.}
}
\description{
The summary method for local polynomial conditional density bandwidth selection objects.
}
\examples{
n=100
x_data = as.matrix(rnorm(n, mean=0, sd=1))
y_data = as.matrix(rnorm(n, mean=0, sd=1))
y_grid = stats::quantile(y_data, seq(from=0.1, to=0.9, by=0.1))
# bandwidth selection
y_grid = stats::quantile(y_data, seq(from=0.1, to=0.9, by=0.1))
model2 = lpcde::lpbwcde(y_data=y_data, x_data=x_data, x=0, y_grid = y_grid, bw_type = "mse-rot")
summary(model2)

}
\seealso{
\code{\link{lpbwcde}} for data-driven bandwidth selection.

Supported methods: \code{\link{coef.lpbwcde}}, \code{\link{print.lpbwcde}}, \code{\link{summary.lpbwcde}}.
}
\author{
Matias D. Cattaneo, Princeton University. \email{cattaneo@princeton.edu}.

Rajita Chandak (maintainer), Princeton University. \email{rchandak@princeton.edu}.

Michael Jansson, University of California Berkeley. \email{mjansson@econ.berkeley.edu}.

Xinwei Ma, University of California San Diego. \email{x1ma@ucsd.edu}.
}
