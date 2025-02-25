################################################################################
#' Print method for local polynomial conditional density bandwidth selection
#'
#' @description The print method for local polynomial conditional density bandwidth selection objects.
#'
#' @param x Class "lpbwcde" object, obtained by calling \code{\link{lpbwcde}}.
#' @param ... Other arguments.
#'
#' @return
#' \item{Display output}{A list of specified options provided to the function.}
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{cattaneo@princeton.edu}.
#'
#' Rajita Chandak (maintainer), Princeton University. \email{rchandak@princeton.edu}.
#'
#' Michael Jansson, University of California Berkeley. \email{mjansson@econ.berkeley.edu}.
#'
#' Xinwei Ma, University of California San Diego. \email{x1ma@ucsd.edu}.
#'
#' @seealso \code{\link{lpbwcde}} for data-driven bandwidth selection.
#'
#' Supported methods: \code{\link{coef.lpbwcde}}, \code{\link{print.lpbwcde}}, \code{\link{summary.lpbwcde}}.
#'
#'
#' @examples
#' n <- 100
#' x_data <- as.matrix(rnorm(n, mean = 0, sd = 1))
#' y_data <- as.matrix(rnorm(n, mean = 0, sd = 1))
#' y_grid <- stats::quantile(y_data, seq(from = 0.1, to = 0.9, by = 0.1))
#' # bandwidth selection
#' y_grid <- stats::quantile(y_data, seq(from = 0.1, to = 0.9, by = 0.1))
#' model2 <- lpcde::lpbwcde(y_data = y_data, x_data = x_data, x = 0,
#'                          y_grid = y_grid, bw_type = "mse-rot")
#' print(model2)
#'
#' @export
print.lpbwcde <- function(x, ...) {
  cat("Call: lpbwcde\n\n")

  cat(paste("Sample size                                      ", x$opt$n, "\n", sep = ""))
  cat(paste("Number of grid points                            ", x$opt$ng, "\n", sep = ""))
  cat(paste("Polynomial order for Y point estimation  (p=)    ", x$opt$p, "\n", sep = ""))
  cat(paste("Polynomial order for X point estimation  (q=)    ", x$opt$q, "\n", sep = ""))
  cat(paste("Order of derivative estimated for Y     (mu=)    ", x$opt$mu, "\n", sep = ""))
  cat(paste("Order of derivative estimated for X     (nu=)    ", x$opt$nu, "\n", sep = ""))
  cat(paste("Kernel function                                  ", x$opt$kernel, "\n", sep = ""))
  cat(paste("Bandwidth method                                 ", x$opt$bw_type, "\n", sep = ""))
  cat("\n")

  cat("Use summary(...) to show bandwidths.\n")
}

################################################################################
#' Summary method for local polynomial conditional density bandwidth selection
#'
#' @description The summary method for local polynomial conditional density bandwidth selection objects.
#'
#' @param object Class "lpbwcde" object, obtained by calling \code{\link{lpbwcde}}.
#' @param ... Additional options, including (i) \code{y_grid} specifies a subset of y_grid points
#'   to display the bandwidth; (ii) \code{gridIndex} specifies the indices of y_grid points
#'   to display the bandwidth.
#'
#' @return
#' \item{Display output}{A list of specified options and a matrix of grid points, bandwidth, and effective sample size.}
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{cattaneo@princeton.edu}.
#'
#' Rajita Chandak (maintainer), Princeton University. \email{rchandak@princeton.edu}.
#'
#' Michael Jansson, University of California Berkeley. \email{mjansson@econ.berkeley.edu}.
#'
#' Xinwei Ma, University of California San Diego. \email{x1ma@ucsd.edu}.
#'
#' @seealso \code{\link{lpbwcde}} for data-driven bandwidth selection.
#'
#' Supported methods: \code{\link{coef.lpbwcde}}, \code{\link{print.lpbwcde}}, \code{\link{summary.lpbwcde}}.
#'
#' @examples
#' n <- 100
#' x_data <- as.matrix(rnorm(n, mean = 0, sd = 1))
#' y_data <- as.matrix(rnorm(n, mean = 0, sd = 1))
#' y_grid <- stats::quantile(y_data, seq(from = 0.1, to = 0.9, by = 0.1))
#' # bandwidth selection
#' y_grid <- stats::quantile(y_data, seq(from = 0.1, to = 0.9, by = 0.1))
#' model2 <- lpcde::lpbwcde(y_data = y_data, x_data = x_data, x = 0,
#'                          y_grid = y_grid, bw_type = "mse-rot")
#' summary(model2)
#'
#' @export
summary.lpbwcde <- function(object, ...) {
  x <- object
  args <- list(...)

  if (is.null(args[["sep"]])) {
    sep <- 5
  } else {
    sep <- args[["sep"]]
  }

  if (is.null(args[["y_grid"]]) & is.null(args[["gridIndex"]])) {
    gridIndex <- 1:nrow(x$BW)
  } else if (is.null(args[["y_grid"]]) & !is.null(args[["gridIndex"]])) {
    gridIndex <- args[["gridIndex"]]
    if (is.null(gridIndex)) {
      gridIndex <- 1:nrow(x$BW)
    } else if (!all(gridIndex %in% 1:nrow(x$BW))) {
      stop(paste("Option gridIndex incorrectly specified. Should be integers between 1 and ", nrow(x$BW), ".\n", sep = ""))
    }
  } else {
    y_grid <- args[["y_grid"]]
    if (is.null(y_grid)) {
      gridIndex <- 1:nrow(x$BW)
    } else if (!is.numeric(y_grid)) {
      stop("Option y_grid incorrectly specified.\n")
    } else {
      gridIndex <- rep(NA, length(y_grid))
      if (min(y_grid) < min(x$BW[, "y_grid"]) | max(y_grid) > max(x$BW[, "y_grid"])) {
        warning("The reporting range exceeds the original estimation range. Option summary(..., y_grid=) should be within the estimation range specified by lpbwcde(..., y_grid=).\n")
      }
      for (j in 1:length(y_grid)) {
        gridIndex[j] <- which.min(abs(x$BW[, "y_grid"] - y_grid[j]))
      }
    }
  }

  cat("Call: lpbwcde\n\n")

  cat(paste("Sample size                                           ", x$opt$n, "\n", sep = ""))
  cat(paste("Polynomial order for Y point estimation      (p=)     ", x$opt$p, "\n", sep = ""))
  cat(paste("Polynomial order for X point estimation      (q=)     ", x$opt$q, "\n", sep = ""))
  if (x$opt$mu == 0) {
    cat(paste("Distribution function estimated             (mu=)   ", x$opt$mu, "\n", sep = ""))
  } else if (x$opt$mu == 1) {
    cat(paste("Density function estimated                   (mu=)    ", x$opt$mu, "\n", sep = ""))
  } else {
    cat(paste("Order of derivative estimated              (mu=)    ", x$opt$mu, "\n", sep = ""))
  }
  cat(paste("Order of derivative estimated for covariates (nu=)    ", x$opt$nu, "\n", sep = ""))
  cat(paste("Kernel function                                       ", x$opt$kernel, "\n", sep = ""))
  cat(paste("Bandwidth method                                      ", x$opt$bw_type, "\n", sep = ""))
  cat("\n")

  ### print output

  cat(paste(rep("=", 12 + 10 + 12), collapse = ""))
  cat("\n")

  cat(format("Index     y_grid", width = 12, justify = "right"))
  cat(format("B.W.", width = 10, justify = "right"))
  cat(format("Eff.n", width = 8, justify = "right"))
  cat("\n")

  cat(paste(rep("=", 12 + 10 + 12), collapse = ""))
  cat("\n")

  jj <- 1
  for (j in gridIndex) {
    cat(format(toString(j), width = 4))
    cat(format(sprintf("%6.4f", x$BW[j, "y_grid"]), width = 10, justify = "right"))
    cat(format(sprintf("%6.4f", x$BW[j, "bw"]), width = 10, justify = "right"))
    cat(format(sprintf("%8.0f", x$BW[j, "nh"]), width = 8, justify = "right"))
    cat("\n")
    if (is.numeric(sep)) {
      if (sep > 0) {
        if (jj %% sep == 0) {
          cat(paste(rep("-", 12 + 10 + 12), collapse = ""))
          cat("\n")
        }
      }
    }
    jj <- jj + 1
  }

  cat(paste(rep("=", 12 + 10 + 12), collapse = ""))
  cat("\n")
}

################################################################################
#' Coef method for local polynomial density bandwidth selection
#'
#' @description The coef method for local polynomial density bandwidth selection objects.
#'
#' @param object Class "lpbwcde" object, obtained by calling \code{\link{lpbwcde}}.
#' @param ... Other arguments.
#'
#' @return
#' \item{Matrix}{A matrix containing y_grid points and selected bandwidths.}
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{cattaneo@princeton.edu}.
#'
#' Rajita Chandak (maintainer), Princeton University. \email{rchandak@princeton.edu}.
#'
#' Michael Jansson, University of California Berkeley. \email{mjansson@econ.berkeley.edu}.
#'
#' Xinwei Ma, University of California San Diego. \email{x1ma@ucsd.edu}.
#'
#' @seealso \code{\link{lpbwcde}} for data-driven bandwidth selection.
#'
#' Supported methods: \code{\link{coef.lpbwcde}}, \code{\link{print.lpbwcde}}, \code{\link{summary.lpbwcde}}.
#'
#' n=100
#' x_data = as.matrix(rnorm(n, mean=0, sd=1))
#' y_data = as.matrix(rnorm(n, mean=0, sd=1))
#' y_grid = stats::quantile(y_data, seq(from=0.1, to=0.9, by=0.1))
#' # bandwidth selection
#' y_grid = stats::quantile(y_data, seq(from=0.1, to=0.9, by=0.1))
#' model2 = lpcde::lpbwcde(y_data=y_data, x_data=x_data, x=0, y_grid = y_grid, bw_type = "mse-rot")
#' coef(model2)
#'
#' @export
coef.lpbwcde <- function(object, ...) {
  object$BW[, c("y_grid", "bw")]
}
