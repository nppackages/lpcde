##################################################################################################
#' @title Data-driven bandwidth selection for local polynomial conditional density estimators
#'
#' @description  \code{\link{lpbwcde}} implements the bandwidth selection methods for local
#'   polynomial based conditionaldensity (and derivatives) estimation proposed and studied
#'   in \insertCite{bernoulli}{lpcde}.
#'
#'   Companion command: \code{\link{lpcde}} for estimation and robust bias-corrected inference.
#'
#'   Related \code{Stata} and \code{R} packages useful for nonparametric estimation and inference are
#'   available at \url{https://nppackages.github.io/}.
#'
#' @param x_data Numeric matrix/data frame, the raw data of covariates.
#' @param y_data Numeric matrix/data frame, the raw data of independent.
#' @param y_grid Numeric, specifies the grid of evaluation points. When set to default, grid points
#'   will be chosen as 0.05-0.95 percentiles of the data, with a step size of 0.05.
#' @param x Numeric, specifies the evaluation point in the x-direction. Default is median of the dataset.
#' @param p Nonnegative integer, specifies the order of the local polynomial for \code{Y} used to
#' construct point estimates. (Default is \code{2}.)
#' @param q Nonnegative integer, specifies the order of the local polynomial for \code{X} used to
#' construct point estimates. (Default is \code{1}.)
#' @param mu Nonnegative integer, specifies the derivative with respect to \code{Y} of the
#' distribution function to be estimated. \code{0} for the distribution function,
#' \code{1} (default) for the density funtion, etc.
#' @param nu Nonnegative integer, specifies the derivative with respect to \code{X} of the
#' distribution function to be estimated.
#' @param grid_spacing String, If equal to "quantile" will generate quantile-spaced grid evaluation points, otherwise will generate equally spaced points.
#' @param ng Int, number of grid points to be used in generating bandwidth estimates.
#' @param kernel_type String, specifies the kernel function, should be one of
#' \code{"triangular"}, \code{"uniform"} or \code{"epanechnikov"}.
#' @param bw_type String, specifies the method for data-driven bandwidth selection. This option will be
#'   ignored if \code{bw} is provided. Implementable with \code{"mse-rot"} (default, mean squared error-optimal
#'   bandwidth selected for each grid point)
  # or (2) \code{"mse-irot"} (integrated mse rule-of-thumb bandwidth with Gaussian
  # reference model).
#' @param regularize Boolean (default TRUE). Option to regularize bandwidth selection to have atleast
#' 20+max(p, q)+1 datapoints when evaluating the estimator.
#'
#' @return
#' \item{BW}{A matrix containing (1) \code{y_grid} (grid point), (2) \code{bw} (bandwidth)}
#' \item{opt}{A list containing options passed to the function.}
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
#' @seealso Supported methods: \code{\link{coef.lpbwcde}},
#' \code{\link{print.lpbwcde}}, \code{\link{summary.lpbwcde}}.
#'
#' @examples
#' # Generate a random sample
#' set.seed(42);
#' x_data = rnorm(2000)
#' y_data = rnorm(2000, mean=x_data)
#' x = 0
#'
#' # Construct bandwidth
#' bw1 <- lpbwcde(y_data = y_data, x_data = x_data, x=x, bw_type = "mse-rot")
#' summary(bw1)
#'
#' # Display bandwidths for a subset of y_grid points
#' summary(bw1, y_grid=bw1$BW[2:5, "y_grid"])
#'
#' @references
#' \insertRef{bernoulli}{lpcde}
#'
#' @export
#'
lpbwcde <- function(y_data, x_data, x, y_grid=NULL, p=NULL, q=NULL, grid_spacing="", ng=NULL,
                    mu=NULL, nu=NULL, kernel_type=c("epanechnikov", "triangular", "uniform"),
                    bw_type=c("mse-rot", "imse-rot"), regularize=NULL){
  ################################################################################
  # Error Checking
  ################################################################################
  # data
  y_data = as.matrix(y_data)
  if (any(is.na(y_data))) {
    warning(paste(sum(is.na(y_data)), " missing ", switch((sum(is.na(y_data))>1)+1, "observation is", "observations are"), " ignored.\n", sep=""))
    y_data = y_data[!is.na(y_data)]
  }

  n = length(y_data)
  if (!is.numeric(y_data) | length(y_data)==0) {
    stop("Data should be numeric, and cannot be empty.\n")
  }

  x_data = as.matrix(x_data)
  if (any(is.na(x_data))) {
    warning(paste(sum(is.na(x_data)), " missing ", switch((sum(is.na(x_data))>1)+1, "observation is", "observations are"), " ignored.\n", sep=""))
    x_data = x_data[!is.na(x_data)]
  }

  d = ncol(x_data)
  if (!is.numeric(x_data) | length(x_data)==0) {
    stop("Data should be numeric, and cannot be empty.\n")
  }

  sd_y = stats::sd(y_data)
  sd_x = apply(x_data, 2, stats::sd)
  mx = apply(x_data, 2, mean)
  my = mean(y_data)
  y_data = (y_data)/sd_y
  x_data = x_data/sd_x
  x = (x-mx)/sd_x
  # y_grid and x_grid
  if (length(y_grid) == 0) {
    if(grid_spacing=="quantile"){
      if (length(bw)>=2){
        y_grid = stats::quantile(y_data, seq(from=0.1, to=0.9, length.out = length(bw)))
        ng = length(y_grid)
      } else{
        if (length(ng)==1){
          y_grid = stats::quantile(y_data, seq(from=0.1, to=0.9, length.out = ng))
        } else{
          y_grid = stats::quantile(y_data, seq(from=0.1, to=0.9, length.out = 19))
          ng =19
        }
      }
    }else{
      gmin = stats::quantile(y_data, probs=0.1)
      gmax = stats::quantile(y_data, probs=0.9)
      if (length(ng)==1){
        y_grid = seq(gmin, gmax, length.out=ng)
      } else{
        y_grid = seq(gmin, gmax, length.out=19)
        ng =19
      }
    }
  } else {
    y_grid = as.vector(y_grid)
    ng = length(y_grid)
    if(!is.numeric(y_grid)) {
      stop("Y grid points should be numeric.\n")
    }
  }

  # evaluation point
  if(length(x)==0){
    flag_no_eval = TRUE
    x = as.vector(apply(x_data, 2, stats::median))
  }else {
    flag_no_eval = FALSE
    x = as.vector(x)
    if(!is.numeric(x)) {
      stop("Evaluation point should be numeric.\n")
    }
  }

  # p
  if (length(p) == 0) {
    if (length(mu) == 0){
      p = 2
    }else {
      p = mu+1
    }
  } else if ((length(p) != 1) | !(p[1]%in%0:20)) {
    stop("Polynomial order p incorrectly specified.\n")
  }

  # q
  if (length(q) == 0) {
    if (length(nu) == 0){
      q = 1
    }else {
      q = nu+1
    }
  } else if ((length(q) != 1) | !(q[1]%in%c(0:20))) {
    stop("Polynomial order (for bias correction) q incorrectly specified.\n")
  }

  # mu
  if (length(mu) == 0) {
    flag_no_mu = TRUE
    mu = min(1, p)
  } else if ((length(mu) != 1) | !(mu[1]%in%c(0:20)) | (mu[1]>p)) {
    stop("Derivative order mu incorrectly specified.\n")
  } else {
    flag_no_mu = FALSE
  }

  # nu
  if (length(nu) == 0) {
    flag_no_nu = TRUE
    nu = min(0, q)
  } else if ((length(nu) != 1) | !(nu[1]%in%c(0:20)) | (nu[1]>q)) {
    stop("Derivative order nu incorrectly specified.\n")
  } else {
    flag_no_nu = FALSE
  }

  # bw_type
  if (length(bw_type) == 0){
    bw_type = "imse-rot"
  } else {
   bw_type = tolower(bw_type)
   bw_type = bw_type[1]
   if (!bw_type%in%c("mse-rot", "imse-rot")){
     stop("Incorrect bandwidth selection method specified.\n")
   }
  }

  # kernel_type
  if (length(kernel_type) == 0) {
    flag_no_kernel = TRUE
    kernel_type = "epanechnikov"
  } else {
    kernel_type = tolower(kernel_type)
    kernel_type = kernel_type[1]
    if (!kernel_type%in%c("triangular", "uniform", "epanechnikov")) {
      stop("Kernel function incorrectly specified.\n")
    } else {
      flag_no_kernel = FALSE
    }
  }

  if (length(regularize) == 0){
    regularize = TRUE
  }
  if(!regularize%in%c(TRUE, FALSE)){
    stop("regularize option must be a Boolean.\n")
  }


  ################################################################################
  # Bandwidth Estimation
  ################################################################################

  if(bw_type == "mse-rot"){
    bw = bw_rot(y_data=y_data, x_data=x_data, y_grid=y_grid, x=x, p=p, q=q, mu=mu, nu=nu, kernel_type=kernel_type, regularize=regularize)

  }else if(bw_type == "imse-rot"){
    bw = bw_irot(y_data=y_data, x_data=x_data, y_grid=y_grid, x=x, p=p, q=q, mu=mu, nu=nu, kernel_type=kernel_type, regularize=regularize)

  }else {
    stop("Invalid bandwidth selection method provided.")

  }


  BW = matrix(NA, ncol=3, nrow=ng)
  BW[, 1] = y_grid
  BW[, 2] = as.vector(bw)
  rownames(BW) = 1:ng
  colnames(BW) = c("y_grid", "bw", "nh")


  if (d == 1){
    for (i in 1:ng) {
      x_idx = which(abs(x_data-x)<=BW[i, 2])
      BW[i, 3] = sum(abs(y_data[x_idx] - BW[i, 1]) <= BW[i, 2])
    }
  } else {
    for (i in 1:ng){
      idx = which(rowSums(abs(sweep(x_data, 2, x))<=BW[i, 2])==d)
      BW[i, 3] = sum(abs(y_data[idx] - BW[i, 1]) <= BW[i, 2])
    }
  }

  Result = list(BW=BW,
                 opt=list(x=x, p=p, q=q, mu=mu, nu=nu, kernel_type=kernel_type, n=n, ng=ng,
                          bw_type=bw_type,
                          data_min=min(y_data), data_max=max(y_data),
                          grid_min=min(y_grid), grid_max=max(y_grid)))

  class(Result) <- c("lpbwcde")

  return (Result)
}
