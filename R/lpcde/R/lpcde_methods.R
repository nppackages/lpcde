#######################################################################################
#' Print method for local polynomial conditional density estimation
#'
#' @description The print method for local polynomial conditional density objects.
#'
#' @param x Class "lpcde" object, obtained from calling \code{\link{lpcde}}.
#' @param ... Additional options.
#'
#' @return
#' \item{Display output}{summary of inputs to \code{lpcde}}
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
#' @seealso \code{\link{lpcde}} for local polynomial conditional density estimation.
#' Supported methods: \code{\link{coef.lpcde}}, \code{\link{confint.lpcde}},
#' \code{\link{plot.lpcde}}, \code{\link{print.lpcde}},
#' \code{\link{summary.lpcde}}, \code{\link{vcov.lpcde}}
#'
#' @examples
#' n=100
#' x_data = as.matrix(rnorm(n, mean=0, sd=1))
#' y_data = as.matrix(rnorm(n, mean=0, sd=1))
#' y_grid = stats::quantile(y_data, seq(from=0.1, to=0.9, by=0.1))
#' # density estimation
#' model1 = lpcde::lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, bw=0.5)
#' print(model1)
#'
#' @export
#'
print.lpcde = function(x, ...){
  cat("Call: lpcde\n\n")

  cat(paste("Sample size                                      ", x$opt$n,        "\n", sep=""))
  cat(paste("Number of grid points                            ", x$opt$ng,       "\n", sep=""))
  cat(paste("Polynomial order for Y point estimation  (p=)    ", x$opt$p,        "\n", sep=""))
  cat(paste("Polynomial order for X point estimation  (q=)    ", x$opt$q,        "\n", sep=""))
  cat(paste("Order of derivative estimated for Y     (mu=)    ", x$opt$mu,        "\n", sep=""))
  cat(paste("Order of derivative estimated for X     (nu=)    ", x$opt$nu,        "\n", sep=""))
  cat(paste("Kernel function                                  ", x$opt$kernel,   "\n", sep=""))
  cat(paste("Bandwidth method                                 ", x$opt$bwselect, "\n", sep=""))
  cat("\n")

  cat("Use summary(...) to show estimates.\n")
}

#######################################################################################
#' Summary method for local polynomial density conditional estimation
#'
#' @description The summary method for local polynomial conditional density objects.
#'
#' @param object Class "lpcde" object, obtained from calling \code{\link{lpcde}}.
#' @param ... Additional options, including (i)\code{y_grid} specifies
#' a subset of grid points in y- directions
#'   to display results; (ii) \code{gridIndex} specifies the indices of grid points
#'   to display results; (iii) \code{alpha} specifies the significance level; (iv)
#'   \code{CIuniform} specifies whether displaying pointwise confidence intervals (\code{FALSE}, default) or
#'   the uniform confidence band (\code{TRUE}); (v) \code{CIsimul} specifies the number of simulations used
#'   to construct critical values (default is \code{2000}).
#'
#' @return
#' \item{Display output}{A list of specified options and a matrix of grid points and estimates.}
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{cattaneo@princeton.edu}.
#'
#' Rajita Chandak (maintainer), Princeton University. \email{rchandak@princeton.edu}
#'
#' Michael Jansson, University of California Berkeley. \email{mjansson@econ.berkeley.edu}.
#'
#' Xinwei Ma, University of California San Diego. \email{x1ma@ucsd.edu}.
#'
#' @seealso \code{\link{lpcde}} for local polynomial conditional density estimation.
#' Supported methods: \code{\link{coef.lpcde}}, \code{\link{confint.lpcde}},
#' \code{\link{plot.lpcde}}, \code{\link{print.lpcde}},
#' \code{\link{summary.lpcde}}, \code{\link{vcov.lpcde}}
#'
#' @examples
#' n=100
#' x_data = as.matrix(rnorm(n, mean=0, sd=1))
#' y_data = as.matrix(rnorm(n, mean=0, sd=1))
#' y_grid = stats::quantile(y_data, seq(from=0.1, to=0.9, by=0.1))
#' # density estimation
#' model1 = lpcde::lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, bw=0.5)
#' summary(model1)
#'
#' @export
summary.lpcde = function(object, ...){
  x = object
  args = list(...)
  if (is.null(args[['alpha']])) { alpha = 0.05 } else { alpha = args[['alpha']] }
  if (is.null(args[['sep']]))   { sep = 5 } else { sep = args[['sep']] }
  if (is.null(args[['CIuniform']]))   { CIuniform <- FALSE } else { CIuniform = args[['CIuniform']] }
  if (is.null(args[['CIsimul']]))   { CIsimul <- 2000 } else { sep = args[['CIsimul']] }

  if (is.null(args[['grid']]) & is.null(args[['gridIndex']])) {
    gridIndex = 1:nrow(x$Estimate)
  } else if (is.null(args[['grid']]) & !is.null(args[['gridIndex']])) {
    gridIndex = args[['gridIndex']]
    if (is.null(gridIndex)) {
      gridIndex = 1:nrow(x$Estimate)
    } else if (!all(gridIndex %in% 1:nrow(x$Estimate))) {
      stop(paste("Option gridIndex incorrectly specified. Should be integers between 1 and ", nrow(x$Estimate), ".\n", sep=""))
    }
  } else {
    grid = args[['grid']]
    if (is.null(grid)) {
      gridIndex = 1:nrow(x$Estimate)
    } else if (!is.numeric(grid)) {
      stop("Option grid incorrectly specified.\n")
    } else {
      gridIndex = rep(NA, length(grid))
      if (min(grid) < min(x$Estimate[, "grid"]) | max(grid) > max(x$Estimate[, "grid"])) {
        warning("The reporting range exceeds the original estimation range. Option summary(..., grid=) should be within the estimation range specified by lpcde(..., grid=).\n")
      }
      for (j in 1:length(grid)) {
        gridIndex[j] = which.min(abs(x$Estimate[, "grid"]-grid[j]))
      }
    }
  }

  cat("Call: lpcde\n\n")

  cat(paste("Sample size                                           ", x$opt$n,        "\n", sep=""))
  cat(paste("Polynomial order for Y point estimation      (p=)     ", x$opt$p,        "\n", sep=""))
  cat(paste("Polynomial order for X point estimation      (q=)     ", x$opt$q,        "\n", sep=""))
  if (x$opt$mu == 0) {
    cat(paste("Distribution function estimated             (mu=)   ", x$opt$mu,        "\n", sep=""))
  } else if (x$opt$mu == 1) {
    cat(paste("Density function estimated                   (mu=)    ", x$opt$mu,        "\n", sep=""))
  } else {
    cat(paste("Order of derivative estimated              (mu=)    ", x$opt$mu,        "\n", sep=""))
  }
  cat(paste("Order of derivative estimated for covariates (nu=)    ", x$opt$nu,        "\n", sep=""))
  cat(paste("Kernel function                                       ", x$opt$kernel,   "\n", sep=""))
  cat(paste("Bandwidth method                                      ", x$opt$bwselect, "\n", sep=""))
  cat("\n")

  # compute CI
  if (CIuniform) {
    if (length(CIsimul) == 0) { CIsimul = 2000 }
    if (!is.numeric(CIsimul) | is.na(CIsimul)) {
      warning("Option CIsimul incorrectly specified. Will only plot pointwise confidence intervals.\n")
      CIuniform = FALSE
      z_val = stats::qnorm(1 - alpha/2)
    }else if (ceiling(CIsimul)<2) {
      warning("Option CIsimul incorrectly specified. Will only plot pointwise confidence intervals.\n")
      CIuniform = FALSE
      z_val = stats::qnorm(1 - alpha/2)
    }else {
      CIsimul = ceiling(CIsimul)
      corrMat = sweep(sweep(x$CovMat$CovMat_RBC, MARGIN=1, FUN="*", STATS=1/x$Estimate[, "se_RBC"]), MARGIN=2, FUN="*", STATS=1/x$Estimate[, "se_RBC"])
      normalSimu = try(
        MASS::mvrnorm(n=CIsimul, mu=rep(0,nrow(corrMat)), Sigma=Matrix::nearPD(corrMat)$mat),
        silent=TRUE)
      if (is.character(normalSimu)) {
        print(normalSimu)
        warning("Variance-Covariance is not positive semidefinite. Will only plot pointwise confidence intervals.\n")
        CIuniform = FALSE
        z_val = stats::qnorm(1 - alpha/2)
      } else {
        z_val = stats::quantile(apply(normalSimu, MARGIN=1, FUN=function(x) {max(abs(x))}), 1 - alpha)
      }
    }
  }else {
    z_val = stats::qnorm(1 - alpha/2)
  }

  CI_l <- x$Estimate[, "est_RBC"] - x$Estimate[, "se_RBC"] * z_val
  CI_r <- x$Estimate[, "est_RBC"] + x$Estimate[, "se_RBC"] * z_val

  # flagNotAllLocalSampleSizeSame <- !all(x$Estimate[, "nh"] == x$Estimate[, "nhu"])

  ### print output
  cat(paste(rep("=", 14 + 10 + 8 + 10 + 10 + 25), collapse="")); cat("\n")

  cat(format(" ", width= 14 ))
  cat(format(" ", width= 10 ))
  cat(format(" ", width= 8 ))
  cat(format("Point", width= 10, justify="right"))
  cat(format("Std." , width= 10, justify="right"))
  cat(format("Robust B.C." , width=25, justify="centre"))
  cat("\n")

  cat(format("Index     Grid"            , width=14, justify="right"))
  cat(format("B.W."              , width=10, justify="right"))
  cat(format("Eff.n"           , width=8 , justify="right"))
  # if (flagNotAllLocalSampleSizeSame) {
  #   cat(format("Uniq.n"           , width=8 , justify="right"))
  # }
  cat(format("Est."            , width=10, justify="right"))
  cat(format("Error"           , width=10, justify="right"))
  if (CIuniform) {
    cat(format(paste("[ Unif. ", floor((1-alpha)*100), "%", " C.I. ]", sep="")
               , width=25, justify="centre"))
  } else {
    cat(format(paste("[ ", floor((1-alpha)*100), "%", " C.I. ]", sep="")
               , width=25, justify="centre"))
  }
  cat("\n")
  cat(paste(rep("=", 14 + 10 + 8 + 10 + 10 + 25), collapse="")); cat("\n")

  jj <- 1
  for (j in gridIndex) {
    cat(format(toString(j), width=4))
    cat(format(sprintf("%6.4f", x$Estimate[j, "y_grid"]), width=10, justify="right"))
    cat(format(sprintf("%6.4f", x$Estimate[j, "bw"])  , width=10, justify="right"))
    cat(format(sprintf("%8.0f", x$eff_n[j]), width=8 , justify="right"))
    # cat(format(sprintf("%8.0f", x$Estimate[j, "nhu"]) , width=8 , justify="right"))
    cat(format(sprintf("%6.4f", x$Estimate[j, "est"]) , width=10, justify="right"))
    cat(format(paste(sprintf("%6.4f", x$Estimate[j, "se"]), sep=""), width=10, justify="right"))
    cat(format(paste(sprintf("%7.4f", CI_l[j]), " , ", sep="")  , width=14, justify="right"))
    cat(format(paste(sprintf("%7.4f", CI_r[j]), sep=""), width=11, justify="left"))
    cat("\n")
    if (is.numeric(sep)) if (sep > 0) if (jj %% sep == 0) {
      cat(paste(rep("-", 14 + 10 + 8 + 10 + 10 + 25), collapse="")); cat("\n")
    }
    jj <- jj + 1
  }

  cat(paste(rep("=", 14 + 10 + 8 + 10 + 10 + 25), collapse=""));
  cat("\n")
}

#######################################################################################
#' Coef method for local polynomial density conditional estimation
#'
#' @description The coef method for local polynomial conditional density objects.
#'
#' @param object Class "lpcde" object, obtained by calling \code{\link{lpcde}}.
#' @param ... Additional options.
#'
#' @return
#' \item{outputs}{A matrix containing the estimates}
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{cattaneo@princeton.edu}.
#'
#' Rajita Chandak (maintainer), Princeton University. \email{rchandak@princeton.edu}
#'
#' Michael Jansson, University of California Berkeley. \email{mjansson@econ.berkeley.edu}.
#'
#' Xinwei Ma, University of California San Diego. \email{x1ma@ucsd.edu}.
#'
#' @seealso \code{\link{lpcde}} for local polynomial conditional density estimation.
#'
#' Supported methods: \code{\link{coef.lpcde}}, \code{\link{confint.lpcde}},
#' \code{\link{plot.lpcde}}, \code{\link{print.lpcde}},
#' \code{\link{summary.lpcde}}, \code{\link{vcov.lpcde}}
#'
#'
#' @examples
#' n=100
#' x_data = as.matrix(rnorm(n, mean=0, sd=1))
#' y_data = as.matrix(rnorm(n, mean=0, sd=1))
#' y_grid = stats::quantile(y_data, seq(from=0.1, to=0.9, by=0.1))
#' # density estimation
#' model1 = lpcde::lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, bw=0.5)
#' coef(model1)
#'
#' @export
coef.lpcde = function(object, ...) {
  object$Estimate
}
#######################################################################################
#' Vcov method for local polynomial density conditional estimation
#'
#' @title Variance-Covariance
#' @description The vcov method for local polynomial conditional density objects.
#'
#' @param object Class "lpdensity" object, obtained by calling \code{\link{lpcde}}.
#' @param ... Additional options.
#'
#' @return
#' \item{stdErr}{A matrix containing grid points and standard errors using p- and q-th order local polynomials.}
#' \item{CovMat}{The variance-covariance matrix corresponding to \code{est}.}
#' \item{CovMat_RBC}{The variance-covariance matrix corresponding to \code{est_RBC}.}
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
#' @seealso \code{\link{lpcde}} for local polynomial conditional density estimation.
#'
#' Supported methods: \code{\link{plot.lpcde}}, \code{\link{print.lpcde}},
#' \code{\link{summary.lpcde}},
#'
#'
#' @examples
#' n=100
#' x_data = as.matrix(rnorm(n, mean=0, sd=1))
#' y_data = as.matrix(rnorm(n, mean=0, sd=1))
#' y_grid = stats::quantile(y_data, seq(from=0.1, to=0.9, by=0.1))
#' # density estimation
#' model1 = lpcde::lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, bw=0.5)
#' vcov(model1)
#'
#' @export
vcov.lpcde = function(object, ...) {

  if (class(object)[1] == "lpbwcde") stop("The vcov method does not support \"lpbwcde\" objects.\n")

  return(list(stdErr=object$Estimate[, c("y_grid", "se", "se_RBC")], CovMat=object$CovMat$CovMat, CovMat_RBC=object$CovMat$CovMat_RBC))
}

#######################################################################################
#' Confint method for local polynomial density conditional estimation
#'
#' @description The confint method for local polynomial conditional density objects.
#'
#' @param object Class "lpdensity" object, obtained by calling \code{\link{lpcde}}.
#' @param parm Integer, indicating which parameters are to be given confidence intervals.
#' @param level Numeric scalar between 0 and 1, the significance level for computing
#' confidence intervals
#' @param alpha Numeric scalar between 0 and 1, specifies the significance level for plotting
#'   confidence intervals/bands.
#' @param CIuniform \code{TRUE} or \code{FALSE} (default), plotting either pointwise confidence intervals (\code{FALSE}) or
#'   uniform confidence bands (\code{TRUE}).
#' @param CIsimul Positive integer, specifies the number of simulations used to construct critical values (default is \code{2000}). This
#'   option is ignored if \code{CIuniform=FALSE}.
#' @param ... Additional options, including (i) \code{grid} specifies a subset of grid points
#'   to display the bandwidth; (ii) \code{gridIndex} specifies the indices of grid points
#'   to display the bandwidth (this is the same as \code{parm});(iii)
#'   \code{CIuniform} specifies whether displaying pointwise confidence intervals
#'   (\code{FALSE}, default) or
#'   the uniform confidence band (\code{TRUE}); (iv) \code{CIsimul} specifies the number of
#'   simulations used to construct critical values (default is 2000).
#'
#' @return
#' \item{Estimate}{A matrix containing grid points, estimates and confidence interval end points using p- and q-th order local polynomials
#' as well as bias-corrected estimates and corresponding confidence intervals.}
#' \item{crit_val}{The critical value used in computing the confidence interval end points.}
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
#' @seealso \code{\link{lpcde}} for local polynomial conditional density estimation.
#'
#' Supported methods: \code{\link{coef.lpcde}}, \code{\link{confint.lpcde}},
#' \code{\link{plot.lpcde}}, \code{\link{print.lpcde}},
#' \code{\link{summary.lpcde}}, \code{\link{vcov.lpcde}}
#'
#'
#' @examples
#' n=100
#' x_data = as.matrix(rnorm(n, mean=0, sd=1))
#' y_data = as.matrix(rnorm(n, mean=0, sd=1))
#' y_grid = stats::quantile(y_data, seq(from=0.1, to=0.9, by=0.1))
#' # density estimation
#' model1 = lpcde::lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, bw=0.5)
#' confint(model1)
#'
#' @export
confint.lpcde <- function(object, parm = NULL, level = NULL, CIuniform=FALSE, CIsimul=2000, alpha=0.05, ...){
  x = object
  if (class(x)[1] == "lpbwcde") {
    stop("The confint method does not support \"lpbwcde\" objects.\n")
  }
  args = list(...)

  if (!is.null(parm)) { args[['grid']] = parm }
  if (!is.null(level)) { args[['alpha']] = 1 - level }

  if (is.null(args[['alpha']])) { alpha = 0.05 } else { alpha = args[['alpha']] }
  if (is.null(args[['CIuniform']]))   { CIuniform = FALSE } else { CIuniform = args[['CIuniform']] }
  if (is.null(args[['CIsimul']]))   { CIsimul = 2000 } else { sep = args[['CIsimul']] }

  if (is.null(args[['grid']]) & is.null(args[['gridIndex']])) {
    gridIndex <- 1:nrow(x$Estimate)
  } else if (is.null(args[['grid']]) & !is.null(args[['gridIndex']])) {
    gridIndex <- args[['gridIndex']]
    if (is.null(gridIndex)) {
      gridIndex <- 1:nrow(x$Estimate)
    } else if (!all(gridIndex %in% 1:nrow(x$Estimate))) {
      stop("Option gridIndex incorrectly specified.\n")
    }
  } else {
    grid <- args[['grid']]
    if (is.null(grid)) {
      gridIndex <- 1:nrow(x$Estimate)
    } else if (!is.numeric(grid)) {
      stop("Option grid incorrectly specified.\n")
    } else {
      gridIndex <- rep(NA, length(grid))
      for (j in 1:length(grid)) {
        gridIndex[j] <- which.min(abs(x$Estimate[, "y_grid"]-grid[j]))
      }
    }
  }

  Estimate = matrix(NA, nrow=length(gridIndex), ncol=7)
  Estimate[, 1] = x$Estimate[gridIndex, "y_grid"]
  Estimate[, 2] = x$Estimate[gridIndex, "est"]
  Estimate[, 3] = x$Estimate[gridIndex, "est_RBC"]

  # critical valus calculations
  if (CIuniform) {
    if (length(CIsimul) == 0) { CIsimul = 2000 }
    if (!is.numeric(CIsimul) | is.na(CIsimul)) {
      warning("Option CIsimul incorrectly specified. Will only plot pointwise confidence intervals.\n")
      CIuniform = FALSE
      z_val = stats::qnorm(1 - alpha/2)
    }else if (ceiling(CIsimul)<2) {
      warning("Option CIsimul incorrectly specified. Will only plot pointwise confidence intervals.\n")
      CIuniform = FALSE
      z_val = stats::qnorm(1 - alpha/2)
    }else {
      CIsimul = ceiling(CIsimul)
      corrMat = sweep(sweep(x$CovMat$CovMat_RBC, MARGIN=1, FUN="*", STATS=1/x$Estimate[, "se_RBC"]), MARGIN=2, FUN="*", STATS=1/x$Estimate[, "se_RBC"])
      normalSimu = try(
        MASS::mvrnorm(n=CIsimul, mu=rep(0,nrow(corrMat)), Sigma=Matrix::nearPD(corrMat)$mat),
        silent=TRUE)
      if (is.character(normalSimu)) {
        print(normalSimu)
        warning("Variance-Covariance is not positive semidefinite. Will only plot pointwise confidence intervals.\n")
        CIuniform = FALSE
        z_val = stats::qnorm(1 - alpha/2)
      } else {
        z_val = stats::quantile(apply(normalSimu, MARGIN=1, FUN=function(x) {max(abs(x))}), 1 - alpha)
      }
    }
  }else {
    z_val = stats::qnorm(1 - alpha/2)
  }


  Estimate[, 4] <- x$Estimate[gridIndex, "est"] - z_val * x$Estimate[gridIndex, "se"]
  Estimate[, 5] <- x$Estimate[gridIndex, "est"] + z_val * x$Estimate[gridIndex, "se"]

  Estimate[, 6] <- x$Estimate[gridIndex, "est_RBC"] - z_val * x$Estimate[gridIndex, "se_RBC"]
  Estimate[, 7] <- x$Estimate[gridIndex, "est_RBC"] + z_val * x$Estimate[gridIndex, "se_RBC"]

  colnames(Estimate) <- c("y_grid", "Est", "Est_RBC", "CI_l", "CI_r", "CI_l_RBC", "CI_r_RBC")

  returnlist = list("Estimate" = Estimate, "crit_val" = z_val)

  return(returnlist)

}

#######################################################################################
#' @title Plot method for local polynomial density conditional estimation
#'
#' @description The plot method for local polynomial density objects.
#' A standard \code{ggplot2} object is returned, hence can be used for further customization.
#'
#' @param ... Class "lpcde" object, obtained from calling \code{\link{lpcde}}.
#' @param alpha Numeric scalar between 0 and 1, specifies the significance level for plotting
#'   confidence intervals/bands.
#' @param type String, one of \code{"line"} (default), \code{"points"} and \code{"both"},
#' specifies how the point estimates are plotted. If more than one is provided,
#' they will be applied to each data series accordingly.
#' @param lty Line type for point estimates, only effective if \code{type} is \code{"line"} or
#'   \code{"both"}. \code{1} for solid line, \code{2} for dashed line, \code{3} for dotted line.
#'   For other options, see the instructions for \code{\link{ggplot2}} . If
#'   more than one is provided, they will be applied to each data series accordingly.
#' @param lwd Line width for point estimates, only effective if \code{type} is \code{"line"} or
#'   \code{"both"}. Should be strictly positive. For other options, see the instructions for
#'   \code{\link{ggplot2}} . If more than one is provided, they will be applied
#'   to each data series accordingly.
#' @param lcol Line color for point estimates, only effective if \code{type} is \code{"line"} or
#'   \code{"both"}. \code{1} for black, \code{2} for red, \code{3} for green, \code{4} for blue.
#'   For other options, see the instructions for \code{\link{ggplot2}} . If
#'   more than one is provided, they will be applied to each data series
#'   accordingly.
#' @param pty Scatter plot type for point estimates, only effective if \code{type} is \code{"points"} or
#'   \code{"both"}. For options, see the instructions for \code{\link{ggplot2}} . If
#'   more than one is provided, they will be applied to each data series
#'   accordingly.
#' @param pwd Scatter plot size for point estimates, only effective if \code{type} is \code{"points"} or
#'   \code{"both"}. Should be strictly positive. If more than one is provided, they will be applied to each data series
#'   accordingly.
#' @param pcol Scatter plot color for point estimates, only effective if \code{type} is \code{"points"} or
#'   \code{"both"}. \code{1} for black, \code{2} for red, \code{3}
#'   for green, \code{4} for blue.
#'   For other options, see the instructions for \code{\link{ggplot2}} . If
#'   more than one is provided, they will be applied to each data series
#'   accordingly.
#' @param y_grid Numeric vector, specifies a subset of grid points
#'   to plot point estimates. This option is effective only if \code{type} is \code{"points"} or
#'   \code{"both"}; or if \code{CItype} is \code{"ebar"} or
#'   \code{"all"}.
#' @param CItype String, one of \code{"region"} (shaded region, default), \code{"line"} (dashed lines),
#'   \code{"ebar"} (error bars), \code{"all"} (all of the previous) or \code{"none"} (no confidence region),
#'   how the confidence region should be plotted. If more than one is provided, they will be applied to each data series
#'   accordingly.
#' @param CIuniform \code{TRUE} or \code{FALSE} (default), plotting either pointwise confidence intervals (\code{FALSE}) or
#'   uniform confidence bands (\code{TRUE}).
#' @param rbc \code{TRUE} or \code{FALSE} (default), plotting confidence intervals and bands with
#' standard estimates (\code{FALSE}) or RBC estimates (\code{TRUE}).
#' @param CIsimul Positive integer, specifies the number of simulations used to construct critical values (default is \code{2000}). This
#'   option is ignored if \code{CIuniform=FALSE}.
#' @param CIshade Numeric, specifies the opaqueness of the confidence region, should be between 0 (transparent) and
#'   1. Default is 0.2. If more than one is provided, they will be applied to each data series
#'   accordingly.
#' @param CIcol Color of the confidence region. \code{1} for black, \code{2} for red, \code{3}
#'   for green, \code{4} for blue.
#'   For other options, see the instructions for \code{\link{ggplot2}} . If
#'   more than one is provided, they will be applied to each data series
#'   accordingly.
#' @param title,xlabel,ylabel Strings, specifies the title of the plot and labels for the x- and y-axis.
#' @param legendTitle String, specifies the legend title.
#' @param legendGroups String vector, specifies the group names used in legend.
#'
#' @return
#' \item{Figure}{A standard \code{ggplot2} object is returned, hence can be used for further customization.}
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{cattaneo@princeton.edu}.
#'
#' Rajita Chandak (maintainer), Princeton University. \email{rchandak@princeton.edu}
#'
#' Michael Jansson, University of California Berkeley. \email{mjansson@econ.berkeley.edu}.
#'
#' Xinwei Ma, University of California San Diego. \email{x1ma@ucsd.edu}.
#'
#' @seealso \code{\link{lpcde}} for local polynomial density estimation.
#' Supported methods: \code{\link{coef.lpcde}}, \code{\link{confint.lpcde}},
#' \code{\link{plot.lpcde}}, \code{\link{print.lpcde}},
#' \code{\link{summary.lpcde}}, \code{\link{vcov.lpcde}}
#'
#' @export
plot.lpcde = function(..., alpha=NULL,type=NULL, lty=NULL, lwd=NULL, lcol=NULL,
                      pty=NULL, pwd=NULL, pcol=NULL, y_grid=NULL,CItype=NULL,
                      CIuniform=FALSE, CIsimul=2000, CIshade=NULL, CIcol=NULL,
                      title=NULL, xlabel=NULL, ylabel=NULL,
                      legendTitle=NULL, legendGroups=NULL, rbc=FALSE) {

  ########################################
  # check how many series are passed in
  ########################################

  x = list(...)
  nfig = length(x)
  if (nfig == 0) {
    stop("Nothing to plot.\n")
  }

  flagToPlot = rep(TRUE, nfig)
  # looping through each series
  for (i in 1:length(x)) {
    # check if is a lpbwcde object
    if (class(x[[i]])[1] == "lpbwcde") {
      flagToPlot[i] = FALSE
      warning(paste("Input ", i, " is an \"lpbwcde\" object, which is not supported by the plot method.\n", sep=""))
      next
    }
    # check if there is only one grid point
    if (nrow(x[[i]]$Estimate) < 2) {
      flagToPlot[i] = FALSE
      warning(paste("At least two grid points are needed to plot input ", i, ".\n", sep=""))
      next
    }
  }

  # select on series that can be plotted
  x = x[flagToPlot]
  nfig = length(x)
  if (nfig == 0) {
    stop("Nothing to plot.\n")
    }

  ############################################################################
  # error handling
  ############################################################################
  # alpha
  if (length(alpha) == 0) {
    alpha = rep(0.05, nfig)
  } else if (!all(alpha>0 & alpha<1)) {
    stop("Significance level incorrectly specified.\n")
  } else {
    alpha = rep(alpha, length.out=nfig)
  }

  # plot type
  if (length(type) == 0) {
    type = rep("line", nfig)
  } else {
    if (!all(type%in%c("line", "points", "both"))) {
      stop("Plotting type incorrectly specified.\n")
    }
    type = rep(type, length.out=nfig)
  }

  # CI type
  if (length(CItype) == 0) {
    CItype = rep("region", nfig)
  } else {
    if (!all(CItype%in%c("region", "line", "ebar", "all", "none"))) {
      stop("Confidence interval type incorrectly specified.\n")
    }
    CItype = rep(CItype, length.out=nfig)
  }

  # line style
  if (length(lty) == 0) {
    lty = rep(1, nfig)
  } else {
    lty = rep(lty, length.out=nfig)
  }
  # line width
  if (length(lwd) == 0) {
    lwd = rep(0.5, nfig)
  } else {
    lwd = rep(lwd, length.out=nfig)
  }
  # line color
  if (length(lcol) == 0) {
    lcol = 1:nfig
  } else {
    lcol = rep(lcol, length.out=nfig)
  }

  # point style
  if (length(pty) == 0) {
    pty = rep(1, nfig)
  } else {
    pty = rep(pty, length.out=nfig)
  }
  # point width
  if (length(pwd) == 0) {
    pwd = rep(1, nfig)
  } else {
    pwd = rep(pwd, length.out=nfig)
  }
  # point color
  if (length(pcol) == 0) {
    pcol = lcol
  } else {
    pcol = rep(pcol, length.out=nfig)
  }

  # CI shade
  if (length(CIshade) == 0) {
    CIshade = rep(0.2, nfig)
  } else {
    CIshade = rep(CIshade, length.out=nfig)
  }
  # CI color
  if (length(CIcol) == 0) {
    CIcol = lcol
  } else {
    CIcol = rep(CIcol, length.out=nfig)
  }

  # legend
  if (length(legendTitle) == 0) {
    legendTitle = ""
  } else {
    legendTitle = legendTitle[1]
  }
  # legend Groups
  if (length(legendGroups) > 0) {
    legendGroups = rep(legendGroups, length.out=nfig)
    legend_default = FALSE
  } else {
    legend_default = TRUE
  }

  # y_grid
  if (!is.null(y_grid)) {
    if (!is.numeric(y_grid)) {
      stop("Option grid incorrectly specified.\n")
    }
  }




  # # x_grid
  # if (!is.null(x_grid)) {
  #   if (!is.numeric(x_grid)) {
  #     stop("Option grid incorrectly specified.\n")
  #   }
  # }

  ########################################
  # initializing plot
  ########################################

  temp_plot = ggplot2::ggplot() + ggplot2::theme_bw()

  CI_l = CI_r = est = Sname = NULL
  ########################################
  # looping over input models
  ########################################
  # estimation range
  estRangeL = estRangeR = c()

  col_all = lty_all = pty_all = v_all = c()

  for (i in 1:nfig) {
    # get derivative order
    v_all <- c(v_all, x[[i]]$opt$mu)
    # get ploting indices
    if (is.null(y_grid)) {
      plotIndex <- 1:nrow(x[[i]]$Estimate)
    }else {
      # selecting grid points corresponding to estimate values computed
      gridTemp <- y_grid[y_grid >= min(x[[i]]$Estimate[, "y_grid"]) & y_grid <= max(x[[i]]$Estimate[, "y_grid"])]
      if (length(gridTemp) == 0) {
        plotIndex <- NULL
      } else {
        plotIndex <- rep(NA, length(gridTemp))
        for (ii in 1:length(gridTemp)) {
          # find minimum distance between grid point and grid values in lpcde object
          plotIndex[ii] <- which.min(abs(gridTemp[ii]-x[[i]]$Estimate[, "y_grid"]))
        }
      }
    }

    # initialize lower and upper bound limits
    estRangeL = min(estRangeL, min(x[[i]]$Estimate[, "y_grid"]))
    estRangeR = max(estRangeR, max(x[[i]]$Estimate[, "y_grid"]))

    # dataset to plot
    data_x = data.frame(x[[i]]$Estimate[, c("y_grid", "est", "est_RBC", "se", "se_RBC"), drop=FALSE])

    # critical values calculations
    if (CIuniform) {
      if (length(CIsimul) == 0) {
        CIsimul = 2000
      }
      if (!is.numeric(CIsimul) | is.na(CIsimul)) {
        warning("Option CIsimul incorrectly specified. Will only plot pointwise confidence intervals.\n")
        CIuniform = FALSE
        z_val = stats::qnorm(1 - alpha/2)
      }else if (ceiling(CIsimul)<2) {
        warning("Option CIsimul incorrectly specified. Will only plot pointwise confidence intervals.\n")
        CIuniform = FALSE
        z_val = stats::qnorm(1 - alpha/2)
      }else {
        CIsimul = ceiling(CIsimul)
        corrMat = sweep(sweep(x[[i]]$CovMat$CovMat_RBC, MARGIN=1, FUN="*", STATS=1/x[[i]]$Estimate[, "se_RBC"], check.margin = FALSE), MARGIN=2, FUN="*", STATS=1/x[[i]]$Estimate[, "se_RBC"], check.margin = FALSE)
        normalSimu = try(
          MASS::mvrnorm(n=CIsimul, mu=rep(0,nrow(corrMat)), Sigma=Matrix::nearPD(corrMat)$mat),
          silent=TRUE)
        if (is.character(normalSimu)) {
          print(normalSimu)
          warning("Variance-Covariance is not positive semidefinite. Will only plot pointwise confidence intervals.\n")
          CIuniform = FALSE
          z_val = stats::qnorm(1 - alpha/2)
        } else {
          z_val = stats::quantile(apply(normalSimu, MARGIN=1, FUN=function(x) {max(abs(x))}), 1 - alpha)
        }
        z_val = stats::confint(x[[i]], CIuniform=TRUE)$crit_val
      }
    }else {
      z_val = stats::qnorm(1 - alpha/2)
    }

    # computing and saving lower and upper confidence interval values
    if(rbc){
      # use RBC values for computation
      data_x$CI_l = data_x$est_RBC - z_val * data_x$se_RBC
      data_x$CI_r = data_x$est_RBC + z_val * data_x$se_RBC
    } else {
      # use standard values for computation
      data_x$CI_l = data_x$est - z_val * data_x$se
      data_x$CI_r = data_x$est + z_val * data_x$se
    }
    # valid estimates
    n_suff = 30
    suff_idx = x[[i]]$eff_n<n_suff
    data_x$CI_l[suff_idx] = 0
    data_x$CI_r[suff_idx] = 0

    # adding legend information to dataset
    if (legend_default) {
      data_x$Sname = paste("Series", i, sep=" ")
      legendGroups = c(legendGroups, data_x$Sname)
    } else {
      data_x$Sname = legendGroups[i]
    }

    ########################################
    # add CI regions to the plot
    if (CItype[i]%in%c("region", "all")){
      temp_plot = temp_plot + ggplot2::geom_ribbon(data=data_x, ggplot2::aes(x=y_grid, ymin=CI_l, ymax=CI_r), alpha=CIshade[i], fill=CIcol[i])
    }

    ########################################
    # add CI lines to the plot
    if (CItype[i]%in%c("line", "all")){
      temp_plot = temp_plot + ggplot2::geom_line(data=data_x, ggplot2::aes(x=y_grid, y=CI_l),
                                                  linetype=2, alpha=1, col=CIcol[i]) +
        ggplot2::geom_line(data=data_x, ggplot2::aes(x=y_grid, y=CI_r), linetype=2,
                           alpha=1, col=CIcol[i])
    }

    ########################################
    # add error bars to the plot
    if (CItype[i]%in%c("ebar", "all") & !is.null(plotIndex)){
      temp_plot = temp_plot + ggplot2::geom_errorbar(data=data_x[plotIndex, ],
                                                      ggplot2::aes(x=y_grid, ymin=CI_l, ymax=CI_r),
                                                      alpha=1, col=CIcol[i], linetype=1)
    }

    ########################################
    # add lines to the plot
    if (type[i]%in%c("line", "both")) {
      temp_plot = temp_plot + ggplot2::geom_line(data=data_x,
                                                 ggplot2::aes(x=y_grid, y=est, colour=Sname,
                                                              linetype=Sname), size=lwd[i])
    }

    ########################################
    # add points to the plot
    if (type[i]%in%c("points", "both") & !is.null(plotIndex)) {
      temp_plot = temp_plot + ggplot2::geom_point(data=data_x[plotIndex, ],
                                                  ggplot2::aes(x=y_grid, y=est, colour=Sname,
                                                               shape=Sname), size=pwd[i])
    }

    if (type[i] == "line") {
      col_all <- c(col_all, lcol[i])
      lty_all <- c(lty_all, lty[i])
      pty_all <- c(pty_all, NA)
    } else if (type[i] == "both") {
      col_all <- c(col_all, lcol[i])
      lty_all <- c(lty_all, lty[i])
      pty_all <- c(pty_all, pty[i])
    } else {
      col_all <- c(col_all, pcol[i])
      lty_all <- c(lty_all, NA)
      pty_all <- c(pty_all, pty[i])
    }

  }

  ########################################
  # change color, line type and point shape back, and customize legend
  ########################################
  # New in v0.2.1 to handle legend
  index <- sort.int(legendGroups, index.return=TRUE)$ix
  temp_plot <- temp_plot + ggplot2::scale_color_manual(values = col_all[index]) +
    ggplot2::scale_linetype_manual(values = lty_all[index]) +
    ggplot2::scale_shape_manual(values = pty_all[index]) +
    ggplot2::guides(colour=ggplot2::guide_legend(title=legendTitle)) +
    ggplot2::guides(linetype=ggplot2::guide_legend(title=legendTitle)) +
    ggplot2::guides(shape=ggplot2::guide_legend(title=legendTitle))

  ########################################
  # add title, x and y labs
  ########################################
  if (is.null(ylabel)) {
    if (all(v_all == v_all[1])) {
      if (v_all[1] == 0) {
        ylabel <- "Distribution function"
      } else if (v_all[1] == 1) {
        ylabel <- "Density"
      } else {
        ylabel <- paste("Density derivative (mu=", v_all[1], ")", sep="")
      }
    } else {
      ylabel <- ""
    }
  }

  if (is.null(xlabel)) {
    xlabel <- ""
  }

  if (is.null(title)) {
    title <- ""
  }

  temp_plot <- temp_plot + ggplot2::labs(x=xlabel, y=ylabel) + ggplot2::ggtitle(title) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  # check plotting range vs estimation range
  if (!is.null(y_grid)) {
    if (min(y_grid) < estRangeL | max(y_grid) > estRangeR) {
      warning("The plotting range exceeds the original estimation range. Option plot(..., grid=) should be within the estimation range specified by lpdensity(..., grid=).\n")
    }
  }

  ########################################
  # return the plot
  ########################################
  return (temp_plot)

}
