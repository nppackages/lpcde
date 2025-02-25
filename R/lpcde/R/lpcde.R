############################################################################################
# This file contains code for generating conditional density estimate (External Functions)
############################################################################################

#' @title Local polynomial conditional density estimation
#' @description \code{\link{lpcde}} implements the local polynomial regression based
#' conditional density (and derivatives). The estimator proposed in
#' \insertCite{bernoulli}{lpcde}.
#' Robust bias-corrected inference methods, both pointwise (confidence intervals) and
#' uniform (confidence bands), are also implemented.
#' @param x_data Numeric matrix/data frame, the raw data of covariates.
#' @param y_data Numeric matrix/data frame, the raw data of independent.
#' @param y_grid Numeric, specifies the grid of evaluation points in the y-direction. When set to default, grid
#' points will be chosen as 0.05-0.95 percentiles of the data, with a step size of 0.05 in
#' y-direction.
#' @param x Numeric, specifies the grid of evaluation points in the x-direction. When set to default,
#' the evaluation point will be chosen as the median of the x data. To generate
#' estimates for multiple conditioning values, please loop over the x values and
#' evaluate the lpcde function at each point.
#' @param bw Numeric, specifies the bandwidth used for estimation. Can be (1) a positive
#' scalar (common bandwidth for all grid points); or (2) a positive numeric vector/matrix
#' specifying bandwidths for each grid point (should be the same dimension as \code{grid}).
#' @param p Nonnegative integer, specifies the order of the local polynomial for \code{Y} used to
#' construct point estimates. (Default is \code{2}.)
#' @param q Nonnegative integer, specifies the order of the local polynomial for \code{X} used to
#' construct point estimates. (Default is \code{1}.)
#' @param p_RBC Nonnegative integer, specifies the order of the local polynomial for \code{Y} used to
#' construct bias-corrected point estimates. (Default is \code{p+1}.)
#' @param q_RBC Nonnegative integer, specifies the order of the local polynomial for \code{X} used to
#' construct bias-corrected point estimates. (Default is \code{q+1}.)
#' @param mu Nonnegative integer, specifies the derivative with respect to \code{Y} of the
#' distribution function to be estimated. \code{0} for the distribution function,
#' \code{1} (default) for the density funtion, etc.
#' @param nu Nonnegative integer, specifies the derivative with respect to \code{X} of the
#' distribution function to be estimated. Default value is \code{0}.
#' @param rbc Boolean. TRUE (default) for rbc calcuations, required for valid uniform inference.
#' @param kernel_type String, specifies the kernel function, should be one of
#' \code{"triangular"}, \code{"uniform"}, and \code{"epanechnikov"}(default).
# @param var_type String, specifies the type of variance estimator to be used.
# choose from \code{"ustat"} or \code{"asymp"}.
#' @param bw_type String, specifies the method for data-driven bandwidth selection. This option will be
#'   ignored if \code{bw} is provided. Implementable with \code{"mse-dpi"} (default, mean squared error-optimal
#'   bandwidth selected for each grid point)
# or (2) \code{"mse-rot"} (rule-of-thumb bandwidth with Gaussian
# reference model).
#' @param ng Int, number of grid points to be used. generates evenly space points over the support of the data.
#' @param grid_spacing String, If equal to "quantile" will generate quantile-spaced grid evaluation points, otherwise will generate equally spaced points.
#' @param normalize Boolean, False (default) returns original estimator, True normalizes estimates to integrate to 1.
#' @param nonneg Boolean, False (default) returns original estimator, True returns maximum of estimate and 0.
#' @param cov_flag String, specifies covariance computation. Must be one of "full" (default), "diag" or "off".
#' @return
#' \item{Estimate}{ A matrix containing (1) \code{grid} (grid points),\cr
#' (2) \code{bw} (bandwidths),\cr
#' (3) \code{est} (point estimates with p-th and q-th order local polynomial),\cr
#' (4) \code{est_RBC} (point estimates with p_RBC-th and q_RBC-th order local polynomial),\cr
#' (5) \code{se} (standard error corresponding to \code{est}. Set to NA if cov_flag="off").
#' (6) \code{se_RBC} (standard error corresponding to \code{est_RBC}). Set to NA if cov_flag="off"}
#' \item{CovMat}{The variance-covariance matrix corresponding to \code{est}. Will be 0 if cov_flag="off" or a diagonal matrix if cov_flag="diag".}
#' \item{opt}{A list containing options passed to the function.}
#' @details Bias correction is only used for the construction of confidence intervals/bands, but not for point estimation.
#' The point estimates, denoted by \code{est}, are constructed using local polynomial estimates of order \code{p} and \code{q},
#'  while the centering of the confidence intervals/bands, denoted by \code{est_RBC},
#'  are constructed using local polynomial estimates of order
#'  \code{p_RBC} and \code{q_RBC}. The confidence intervals/bands take the form:
#'  \code{[est_RBC - cv * SE(est_RBC) , est_RBC + cv * SE(est_RBC)]}, where \code{cv} denotes
#'  the appropriate critical value and \code{SE(est_RBC)} denotes an standard error estimate for
#'  the centering of the confidence interval/band. As a result, the confidence intervals/bands
#'  may not be centered at the point estimates because they have been bias-corrected.
#'  Setting \code{p_RBC} equal to \code{p} and \code{q_RBC} to \code{q}, results on centered
#'  at the point estimate confidence intervals/bands, but requires undersmoothing for
#'  valid inference (i.e., (I)MSE-optimal bandwdith for the density point estimator cannot
#'  be used). Hence the bandwidth would need to be specified manually when \code{q=p},
#'  and the point estimates will not be (I)MSE optimal. See Cattaneo, Jansson and Ma
#'  (2020a, 2020b) for details, and also Calonico, Cattaneo, and Farrell (2018, 2020)
#'  for robust bias correction methods.
#'
#'  Sometimes the density point estimates may lie outside
#'  of the confidence intervals/bands, which can happen if the underlying distribution exhibits
#'  high curvature at some evaluation point(s). One possible solution in this case is to
#'  increase the polynomial order \code{p} or to employ a smaller bandwidth.
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
#' @seealso Supported methods: \code{\link{coef.lpcde}}, \code{\link{confint.lpcde}},
#' \code{\link{plot.lpcde}}, \code{\link{print.lpcde}},
#' \code{\link{summary.lpcde}}, \code{\link{vcov.lpcde}}
#'
#'
#' @examples
#' # Density estimation example
#' n <- 500
#' x_data <- matrix(rnorm(n, mean = 0, sd = 1))
#' y_data <- matrix(rnorm(n, mean = x_data, sd = 1))
#' y_grid <- seq(from = -1, to = 1, length.out = 5)
#' model1 <- lpcde::lpcde(x_data = x_data, y_data = y_data, y_grid = y_grid, x = 0, bw = 0.5)
#' # summary of estimation
#' summary(model1)
#'
#' @references
#' \insertRef{bernoulli}{lpcde}\cr
#' \insertRef{JASA}{lpcde}\cr
#' \insertRef{rbc}{lpcde}\cr
#' \insertRef{lpdensitypaper}{lpcde}
#' @export
lpcde <- function(x_data, y_data, y_grid = NULL, x = NULL, bw = NULL, p = NULL, q = NULL,
                  p_RBC = NULL, q_RBC = NULL, mu = NULL, nu = NULL, rbc = TRUE, ng = NULL,
                  cov_flag = c("full", "diag", "off"),
                  normalize = FALSE, nonneg = FALSE, grid_spacing = "",
                  kernel_type = c("epanechnikov", "triangular", "uniform"),
                  bw_type = NULL) {
  ################################################################################
  # Error Checking
  ################################################################################
  # data
  y_data <- as.matrix(y_data)
  if (any(is.na(y_data))) {
    warning(paste(sum(is.na(y_data)), " missing ", switch((sum(is.na(y_data)) > 1) + 1,
      "observation is",
      "observations are"
    ), " ignored.\n", sep = ""))
    y_data <- y_data[!is.na(y_data)]
  }

  n <- length(y_data)
  if (!is.numeric(y_data) | length(y_data) == 0) {
    stop("Data should be numeric, and cannot be empty.\n")
  }

  x_data <- as.matrix(x_data)
  if (any(is.na(x_data))) {
    warning(paste(sum(is.na(x_data)), " missing ", switch((sum(is.na(x_data)) > 1) + 1,
      "observation is",
      "observations are"
    ), " ignored.\n", sep = ""))
    x_data <- x_data[!is.na(x_data)]
  }

  n <- ncol(x_data)
  if (!is.numeric(x_data) | length(x_data) == 0) {
    stop("Data should be numeric, and cannot be empty.\n")
  }

  if (length(cov_flag) == 0) {
    cov_flag <- "full"
  } else {
    cov_flag <- cov_flag[1]
    if (!cov_flag %in% c("full", "diag", "off")) {
      stop("Incorrect covariance estimation flag provided. Please see the documentation on available options.")
    }
  }


  # sd_y = stats::sd(y_data)
  # sd_x = apply(x_data, 2, stats::sd)
  # mx = apply(x_data, 2, mean)
  # my = mean(y_data)
  # y_data = (y_data)/sd_y
  # x_data = x_data/sd_x
  # grid
  if (length(y_grid) == 0) {
    if (grid_spacing == "quantile") {
      if (length(bw) >= 2) {
        y_grid <- stats::quantile(y_data, seq(from = 0.1, to = 0.9, length.out = length(bw)))
        ng <- length(y_grid)
      } else {
        if (length(ng) == 1) {
          y_grid <- stats::quantile(y_data, seq(from = 0.1, to = 0.9, length.out = ng))
        } else {
          y_grid <- stats::quantile(y_data, seq(from = 0.1, to = 0.9, length.out = 19))
          ng <- 19
        }
      }
    } else {
      gmin <- stats::quantile(y_data, probs = 0.1)
      gmax <- stats::quantile(y_data, probs = 0.9)
      if (length(bw) >= 2) {
        y_grid <- seq(gmin, gmax, length.out = length(bw))
        ng <- length(y_grid)
      } else {
        if (length(ng) == 1) {
          y_grid <- seq(gmin, gmax, length.out = ng)
        } else {
          y_grid <- seq(gmin, gmax, length.out = 19)
          ng <- 19
        }
      }
    }
  } else {
    y_grid <- as.vector(y_grid)
    ng <- length(y_grid)
    if (!is.numeric(y_grid)) {
      stop("Y grid points should be numeric.\n")
    }
  }

  # evaluation point
  if (length(x) == 0) {
    x <- as.vector(apply(x_data, 2, stats::median))
  } else {
    x <- as.vector(x)
    if (!is.numeric(x)) {
      stop("Evaluation point should be numeric.\n")
    }
  }

  # p
  if (length(p) == 0) {
    if (length(mu) == 0) {
      p <- 2
    } else {
      p <- mu + 1
    }
  } else if ((length(p) != 1) | !(p[1] %in% 0:20)) {
    stop("Polynomial order p incorrectly specified.\n")
  }

  # q
  if (length(q) == 0) {
    if (length(nu) == 0) {
      q <- 1
    } else {
      q <- nu + 1
    }
  } else if ((length(q) != 1) | !(q[1] %in% c(0:20))) {
    stop("Polynomial order (for bias correction) q incorrectly specified.\n")
  }

  # p_RBC
  if (length(p_RBC) == 0) {
    p_RBC <- p + 1
  } else if ((length(p_RBC) != 1) | !(p_RBC[1] %in% c(0:20)) | (p_RBC[1] < p)) {
    stop("Polynomial order (for bias correction) q incorrectly specified.\n")
  }
  # q_RBC
  if (length(q_RBC) == 0) {
    q_RBC <- q + 1
  } else if ((length(q_RBC) != 1) | !(q_RBC[1] %in% c(0:20)) | (q_RBC[1] < q)) {
    stop("Polynomial order (for bias correction) q incorrectly specified.\n")
  }

  if (p_RBC == p && q_RBC == q) {
    warning("RBC polynomial order same as estimation order. No Bias Correction Implemented.\n")
  }

  # mu
  if (length(mu) == 0) {
    mu <- min(1, p)
  } else if ((length(mu) != 1) | !(mu[1] %in% c(0:20)) | (mu[1] > p)) {
    stop("Derivative order v incorrectly specified.\n")
  }

  # nu
  if (length(nu) == 0) {
    nu <- min(0, q)
  } else if ((length(nu) != 1) | !(nu[1] %in% c(0:20)) | (nu[1] > q)) {
    stop("Derivative order v incorrectly specified.\n")
  }

  # kernel_type
  if (length(kernel_type) == 0) {
    kernel_type <- "epanechnikov"
  } else {
    kernel_type <- tolower(kernel_type)
    kernel_type <- kernel_type[1]
    if (!kernel_type %in% c("triangular", "uniform", "epanechnikov")) {
      stop("Kernel function incorrectly specified.\n")
    }
  }

  # bw
  if (length(bw) == 0) {
    if (length(bw_type) == 0) {
      bw_type <- "imse-rot"
      bw <- lpbwcde(
        y_data = y_data, x_data = x_data, x = x, y_grid = y_grid, p = p, q = q, mu = mu,
        nu = nu, kernel_type = kernel_type, bw_type = bw_type
      )$BW[, 2]
    } else {
      bw_type <- tolower(bw_type[1])
      bw <- lpbwcde(
        y_data = y_data, x_data = x_data, x = x, y_grid = y_grid, p = p, q = q, mu = mu,
        nu = nu, kernel_type = kernel_type, bw_type = bw_type
      )$BW[, 2]
    }
    if (!bw_type %in% c("mse-rot", "imse-rot")) {
      stop("Incorrect bandwidth selection method specified.\n")
    }
  } else if (length(bw) == 1) {
    if (!is.numeric(bw) | bw <= 0) {
      stop("Bandwidth incorrectly specified.\n")
    } else {
      bw <- rep(bw, ng)
      bw_type <- "user provided"
    }
  } else {
    bw <- as.vector(bw)
    if (!is.numeric(bw)) {
      stop("Bandwidth incorrectly specified.\n")
    } else if (length(bw) != ng) {
      stop("Bandwidth has to be the same length as grid.\n")
    } else {
      bw_type <- "user provided"
    }
  }

  ################################################################################
  # Point Estimation and Standard Error
  ################################################################################

  lpcdest <- lpcde_fn(
    y_data = y_data, x_data = x_data, y_grid = y_grid,
    x = x, p = p, q = q, p_RBC = p_RBC, q_RBC = q_RBC,
    bw = bw, mu = mu, nu = nu, cov_flag = cov_flag,
    kernel_type = kernel_type, rbc = rbc
  )
  rownames(lpcdest$est) <- 1:ng

  ################################################################################
  # Normalizing
  ################################################################################

  if (nonneg == TRUE) {
    lpcdest$est[, 3] <- replace(lpcdest$est[, 3], lpcdest$est[, 3] < 0, 0)
  }
  if (normalize == TRUE) {
    grid_diff <- c(diff(y_grid), diff(utils::tail(y_grid, 2)))
    c <- sum(lpcdest$est[, 3] * grid_diff)
    lpcdest$est[, 3] <- lpcdest$est[, 3] / c
  }

  ################################################################################
  # Return
  ################################################################################
  Result <- list(
    Estimate = lpcdest$est,
    CovMat = lpcdest$CovMat,
    eff_n = lpcdest$eff_n,
    opt = list(
      p = p, q = q, p_RBC = p_RBC, q_RBC = q_RBC,
      mu = mu, nu = nu, kernel = kernel_type, n = length(y_data), ng = ng,
      bw_type = bw_type, bw = bw, xeval = x, cov_flag = cov_flag,
      y_data_min = min(y_data), y_data_max = max(y_data),
      x_data_min = min(x_data), x_data_max = max(x_data),
      grid_min = min(y_grid), grid_max = max(y_grid)
    )
  )

  if (any(lpcdest$eff_n <= 5)) {
    warning("Some evaluation points do not have enough data to produce reliable results.")
  }
  if (lpcdest$singular_flag == TRUE) {
    warning("Singular matrices encountered. May affect estimates.")
  }

  class(Result) <- c("lpcde")
  return(Result)
}
