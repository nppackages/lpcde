#######################################################################################
# This file contains code for generating point estimates (Internal Functions)
#######################################################################################

#' @title lpcde_fn: Conditional density estimator.
#' @description Function for estimating the density function and its derivatives.
#' @param y_data Response variable dataset, vector.
#' @param x_data Covariate dataset, vector or matrix.
#' @param y_grid Numeric vector, specifies the grid of evaluation points along y-direction.
#' @param x Numeric vector or matrix, specifies the grid of evaluation points along x-direction.
#' @param p Polynomial order for y.
#' @param q Polynomial order for covariates.
#' @param p_RBC Nonnegative integer, specifies the order of the local polynomial for \code{Y} used to
#' construct bias-corrected point estimates. (Default is \code{p+1}.)
#' @param q_RBC Nonnegative integer, specifies the order of the local polynomial for \code{X} used to
#' construct bias-corrected point estimates. (Default is \code{q+1}.)
#' @param bw Numeric, bandwidth vector.
#' @param mu Degree of derivative with respect to y.
#' @param nu Degree of derivative with respect to x.
#' @param kernel_type Kernel function choice.
#' @param cov_flag Flag for covariance estimation option.
#' @param rbc Boolean for whether to return RBC estimate and standard errors.
# @param var_type String. type of variance estimator to implement.
# choose from "ustat" and "asymp".
#' @return Conditional density estimate at all grid points.
#' @keywords internal
#' @importFrom Rdpack reprompt
lpcde_fn <- function(y_data, x_data, y_grid, x, p, q, p_RBC, q_RBC, bw, mu, nu,
                     cov_flag, kernel_type, rbc = FALSE) {
  sd_y <- stats::sd(y_data)
  sd_x <- apply(x_data, 2, stats::sd)
  mx <- apply(x_data, 2, mean)
  my <- mean(y_data)
  # y_data = (y_data - my)/sd_y
  # y_grid = (y_grid-my)/sd_y
  d <- ncol(x_data)
  x_data <- sweep(x_data, 2, mx) / sd_x
  x <- matrix(x, ncol = d)
  x <- sweep(x, MARGIN = 2, STATS = matrix(mx, ncol = d)) / sd_x
  # initializing output vectors
  est <- matrix(0L, nrow = length(y_grid), ncol = 1)
  se <- matrix(0L, nrow = length(y_grid), ncol = 1)

  # estimate
  f_hat_val <- fhat(
    x_data = x_data, y_data = y_data, x = x, y_grid = y_grid, p = p, q = q,
    mu = mu, nu = nu, h = bw, kernel_type = kernel_type
  )
  est <- f_hat_val$est
  eff.n <- f_hat_val$eff.n
  est_flag <- f_hat_val$singular_flag

  # standard errors
  if (cov_flag == "off") {
    covMat <- NA
    c_flag <- FALSE
    se <- rep(NA, length(y_grid))
  } else {
    covmat <- cov_hat(
      x_data = x_data, y_data = y_data, x = x, y_grid = y_grid, p = p, q = q,
      mu = mu, nu = nu, h = bw, kernel_type = kernel_type, cov_flag = cov_flag
    )
    covMat <- covmat$cov
    c_flag <- covmat$singular_flag
    se <- sqrt(abs(diag(covMat))) * sd_y * mean(sd_x)
  }

  if (rbc) {
    est_rbc <- matrix(0L, nrow = length(y_grid), ncol = 1)
    se_rbc <- matrix(0L, nrow = length(y_grid), ncol = 1)
    if (p_RBC == p && q_RBC == q) {
      est_rbc <- est
      se_rbc <- se
      covMat_rbc <- covMat

      rbc_flag <- f_hat_val$singular_flag
      c_rbc_flag <- covmat$singular_flag
    } else {
      # estimate
      f_hat_rbc <- fhat(
        x_data = x_data, y_data = y_data, x = x, y_grid = y_grid, p = p_RBC,
        q = q_RBC, mu = mu, nu = nu, h = bw, kernel_type = kernel_type
      )

      est_rbc <- f_hat_rbc$est
      rbc_flag <- f_hat_rbc$singular_flag
      # covariance matrix

      # standard errors
      if (cov_flag == "off") {
        covMat_rbc <- NA
        c_rbc_flag <- FALSE
        se_rbc <- rep(NA, length(y_grid))
      } else {
        covmat_rbc <- cov_hat(
          x_data = x_data, y_data = y_data, x = x, y_grid = y_grid,
          p = p_RBC, q = q_RBC, mu = mu, nu = nu, h = bw,
          kernel_type = kernel_type, cov_flag = cov_flag
        )
        covMat_rbc <- covmat_rbc$cov
        c_rbc_flag <- covmat_rbc$singular_flag
        se_rbc <- sqrt(abs(diag(covMat_rbc))) * sd_y * mean(sd_x)
      }
    }

    # with rbc results
    if (est_flag == TRUE | c_flag == TRUE | rbc_flag == TRUE | c_rbc_flag == TRUE) {
      singular_flag <- TRUE
    } else {
      singular_flag <- FALSE
    }
    x <- x * sd_x + mx
    # y_grid = y_grid*sd_y + my
    estimate <- cbind(y_grid, bw, est, est_rbc, se, se_rbc)
    colnames(estimate) <- c("y_grid", "bw", "est", "est_RBC", "se", "se_RBC")
    rownames(estimate) <- c()
    est_result <- list(
      "est" = estimate,
      "CovMat" = list(
        "CovMat" = covMat,
        "CovMat_RBC" = covMat_rbc
      ),
      "x" = x, "eff_n" = eff.n, "singular_flag" = singular_flag
    )
  } else {
    # without rbc results
    # setting rbc values to be the same as non-rbc values
    est_rbc <- est
    se_rbc <- se
    covMat_rbc <- covMat
    # generating matrices and list to return
    x <- x * sd_x + mx
    # y_grid = y_grid*sd_y + my
    estimate <- cbind(y_grid, bw, est, est_rbc, se, se_rbc)
    colnames(estimate) <- c("y_grid", "bw", "est", "est_RBC", "se", "se_RBC")
    rownames(estimate) <- c()
    est_result <- list(
      "est" = estimate,
      "CovMat" = list(
        "CovMat" = covMat,
        "CovMat_RBC" = covMat_rbc
      ),
      "x" = x, "eff_n" = eff.n, "singular_flag" = singular_flag
    )
  }

  return(est_result)
}

#' @title Estimator construction
#' @description Function for estimating the density function and its derivatives.
#' @param y_data Response variable dataset, vector.
#' @param x_data Covariate dataset, vector or matrix.
#' @param y_grid Numeric vector, specifies the grid of evaluation points along y-direction.
#' @param x Numeric vector or matrix, specifies the grid of evaluation points along x-direction.
#' @param p Polynomial order for y.
#' @param q Polynomial order for covariates.
#' @param h Numeric, bandwidth vector.
#' @param mu Degree of derivative with respect to y.
#' @param nu Degree of derivative with respect to x.
#' @param kernel_type Kernel function choice.
#' @return Conditional density estimate at all grid points.
#' @keywords internal
fhat <- function(x_data, y_data, x, y_grid, p, q, mu, nu, h, kernel_type) {
  # setting constants
  n <- length(y_data)
  d <- ncol(x_data)
  ng <- length(y_grid)
  x <- matrix(x, ncol = d)
  singular_flag <- FALSE

  # x basis vector
  e_nu <- basis_vec(x, q, nu)

  # y basis vector
  e_mu <- basis_vec(0, p, mu)

  f_hat <- matrix(0L, nrow = ng)
  nh_vec <- matrix(0L, nrow = ng)

  if (length(unique(h)) == 1) {
    h <- h[1]
    # localization for x
    idx <- which(rowSums(abs(sweep(x_data, 2, x)) <= h) == d)
    # print(x_data)

    x_idx <- matrix(x_data[idx, ], ncol = d)
    y_idx <- y_data[idx]

    # idx of ordering wrt y
    sort_idx <- sort(y_idx, index.return = TRUE)$ix

    # sorting datasets
    y_sorted <- as.matrix(y_idx[sort_idx])
    x_sorted <- matrix(x_idx[sort_idx, ], ncol = d)

    # x constants
    x_scaled <- sweep(x_sorted, 2, x) / (h^d)
    if (length(x_scaled) == 0) {
      bx <- 0
    } else {
      if (check_inv(S_x(x_scaled, q, kernel_type) / (n * h^d))[1] == TRUE) {
        sx_mat <- solve(S_x(x_scaled, q, kernel_type) / (n * h^d))
      } else {
        singular_flag <- TRUE
        sx_mat <- matrix(0L, nrow = length(e_nu), ncol = length(e_nu))
      }
      bx <- b_x(x_scaled, sx_mat, e_nu, q, kernel_type)
    }

    for (j in 1:ng) {
      y <- y_grid[j]
      y_scaled <- (y_sorted - y) / h
      y_elems <- which(abs(y_scaled) <= 1)

      if (length(y_elems) <= 5) {
        f_hat[j] <- 0
        nh_vec[j] <- length(y_elems)
      } else {
        if (mu == 0) {
          if (check_inv(S_x(as.matrix(x_scaled[y_elems, ]), q, kernel_type) / (n * h^d))[1] == TRUE) {
            sx_mat <- solve(S_x(as.matrix(x_scaled[y_elems, ]), q, kernel_type) / (n * h^d))
          } else {
            singular_flag <- TRUE
            sx_mat <- matrix(0L, nrow = length(e_nu), ncol = length(e_nu))
          }
          bx <- b_x(x_scaled, sx_mat, e_nu, q, kernel_type)
          # sy matrix
          if (check_inv(S_x(as.matrix(y_scaled), p, kernel_type) / (n * h))[1] == TRUE) {
            sy_mat <- solve(S_x(as.matrix(y_scaled), p, kernel_type) / (n * h))
          } else {
            singular_flag <- TRUE
            sy_mat <- matrix(0L, nrow = length(e_mu), ncol = length(e_mu))
          }

          # y constants
          ax <- b_x(y_scaled, sy_mat, e_mu, p, kernel_type)

          # adding and multiplying
          f_hat[j] <- ax[y_elems] %*% cumsum(bx[y_elems])

          # number of datapoints used
          nh_vec[j] <- length(y_elems)
        } else {
          # sy matrix
          if (check_inv(S_x(as.matrix(y_scaled[y_elems]), p, kernel_type) / (n * h))[1] == TRUE) {
            sy_mat <- solve(S_x(as.matrix(y_scaled[y_elems]), p, kernel_type) / (n * h))
          } else {
            singular_flag <- TRUE
            sy_mat <- matrix(0L, nrow = length(e_mu), ncol = length(e_mu))
          }

          # y constants
          ax <- b_x(as.matrix(y_scaled[y_elems]), sy_mat, e_mu, p, kernel_type)

          # adding and multiplying
          f_hat[j] <- ax %*% cumsum(bx[y_elems])

          # number of datapoints used
          nh_vec[j] <- length(y_elems)
        }
      }
    }

    # scaling
    f_hat <- f_hat / (n^2 * h^(d + mu + nu + 1))
  } else {
    for (j in 1:ng) {
      # localization for x
      idx <- which(rowSums(abs(sweep(x_data, 2, x)) <= h[j]) == d)

      x_data_loc <- matrix(x_data[idx, ], ncol = d)
      y_data_loc <- y_data[idx]

      # idx of ordering wrt y
      sort_idx <- sort(y_data_loc, index.return = TRUE)$ix

      # sorting datasets
      y_data_loc <- as.matrix(y_data_loc[sort_idx])
      x_data_loc <- matrix(x_data_loc[sort_idx, ], ncol = d)

      # x constants
      x_scaled <- sweep(x_data_loc, 2, x) / (h[j]^d)
      if (check_inv(S_x(as.matrix(x_scaled), q, kernel_type) / (n * h[j]^d))[1] == TRUE) {
        sx_mat <- solve(S_x(as.matrix(x_scaled), q, kernel_type) / (n * h[j]^d))
      } else {
        singular_flag <- TRUE
        sx_mat <- matrix(0L, nrow = length(e_nu), ncol = length(e_nu))
      }
      bx <- b_x(as.matrix(x_scaled), sx_mat, e_nu, q, kernel_type)

      y <- y_grid[j]
      y_scaled <- (y_data_loc - y) / h[j]
      y_elems <- which(abs(y_scaled) <= 1)

      if (length(y_elems) <= 5) {
        f_hat[j] <- 0
        nh_vec[j] <- length(y_elems)
      } else {
        if (mu == 0) {
          if (check_inv(S_x(as.matrix(x_scaled[y_elems, ]), q, kernel_type) / (n * h[j]^d))[1] == TRUE) {
            sx_mat <- solve(S_x(as.matrix(x_scaled[y_elems, ]), q, kernel_type) / (n * h[j]^d))
          } else {
            singular_flag <- TRUE
            sx_mat <- matrix(0L, nrow = length(e_nu), ncol = length(e_nu))
          }
          bx <- b_x(as.matrix(x_scaled), sx_mat, e_nu, q, kernel_type)
          # sy matrix
          if (check_inv(S_x(as.matrix(y_scaled), p, kernel_type) / (n * h[j]))[1] == TRUE) {
            sy_mat <- solve(S_x(as.matrix(y_scaled), p, kernel_type) / (n * h[j]))
          } else {
            singular_flag <- TRUE
            sy_mat <- matrix(0L, nrow = length(e_mu), ncol = length(e_mu))
          }

          # y constants
          ax <- b_x(matrix(y_scaled), sy_mat, e_mu, p, kernel_type)

          # adding and multiplying
          f_hat[j] <- ax %*% cumsum(bx[y_elems])

          # scaling
          f_hat[j] <- f_hat[j] / (n^2 * h[j]^(d + mu + nu + 1))
          # number of datapoints used
          nh_vec[j] <- length(y_elems)
        } else {
          # sy matrix
          if (check_inv(S_x(as.matrix(y_scaled[y_elems]), p, kernel_type) / (n * h[j]))[1] == TRUE) {
            sy_mat <- solve(S_x(as.matrix(y_scaled[y_elems]), p, kernel_type) / (n * h[j]))
          } else {
            singular_flag <- TRUE
            sy_mat <- matrix(0L, nrow = length(e_mu), ncol = length(e_mu))
          }

          # y constants
          ax <- b_x(matrix(y_scaled[y_elems]), sy_mat, e_mu, p, kernel_type)

          # adding and multiplying
          f_hat[j] <- ax %*% cumsum(bx[y_elems])

          # scaling
          f_hat[j] <- f_hat[j] / (n^2 * h[j]^(d + mu + nu + 1))
          # number of datapoints used
          nh_vec[j] <- length(y_elems)
        }
      }
    }
  }
  return(list("est" = f_hat, "eff.n" = nh_vec, "singular_flag" = singular_flag))
}

#' @title cov_hat: covariance estimator
#' @description Function for estimating the variance-covariance matrix.
#' @param y_data Response variable dataset, vector.
#' @param x_data Covariate dataset, vector or matrix.
#' @param y_grid Numeric vector, specifies the grid of evaluation points along
#' y-direction.
#' @param x Numeric vector or matrix, specifies the grid of evaluation points
#' along x-direction.
#' @param p Polynomial order for y.
#' @param q Polynomial order for covariates.
#' @param h Numeric, bandwidth vector.
#' @param mu Degree of derivative with respect to y.
#' @param nu Degree of derivative with respect to x.
#' @param kernel_type Kernel function choice.
#' @param cov_flag Flag for covariance estimation option.
#' @return Covariance matrix for all the grid points
#' @keywords internal
cov_hat <- function(x_data, y_data, x, y_grid, p, q, mu, nu, h, kernel_type, cov_flag) {
  # setting constants
  n <- length(y_data)
  d <- ncol(x_data)
  ng <- length(y_grid)

  singular_flag <- FALSE
  # x basis vector
  e_nu <- basis_vec(x, q, nu)

  # y basis vector
  e_mu <- basis_vec(1, p, mu)

  if (length(unique(h)) == 1) {
    h <- h[1]
    # localization for x
    idx <- which(rowSums(abs(sweep(x_data, 2, x)) <= h) == d)

    x_idx <- as.matrix(x_data[idx, ])
    y_idx <- y_data[idx]

    # idx of ordering wrt y
    sort_idx <- sort(y_idx, index.return = TRUE)$ix

    # sorting datasets
    y_sorted <- as.matrix(y_idx[sort_idx])
    x_sorted <- as.matrix(x_idx[sort_idx, ])

    # x constants
    x_scaled <- sweep(x_sorted, 2, x) / (h^d)
    if (length(x_scaled) == 0) {
      bx <- 0
    } else {
      if (check_inv(S_x(x_scaled, q, kernel_type) / (n * h^d))[1] == TRUE) {
        sx_mat <- solve(S_x(x_scaled, q, kernel_type) / (n * h^d))
      } else {
        singular_flag <- TRUE
        sx_mat <- matrix(0L, nrow = length(e_nu), ncol = length(e_nu))
      }
      bx <- b_x(x_scaled, sx_mat, e_nu, q, kernel_type)
    }


    # initialize matrix
    c_hat <- matrix(0L, nrow = ng, ncol = ng)

    if (cov_flag == "diag") {
      c_hat <- matrix(NA, nrow = ng, ncol = ng)
      for (i in 1:ng) {
        # relevant entries w.r.t. y and y_prime
        y <- y_grid[i]
        y_prime <- y_grid[i]
        y_scaled <- (y_sorted - y) / h
        yp_scaled <- (y_sorted - y_prime) / h
        y_elems <- which(abs(y_scaled) <= 1)
        yp_elems <- which(abs(yp_scaled) <= 1)
        elems <- intersect(y_elems, yp_elems)

        if (length(elems) <= 5) {
          c_hat[i, i] <- 0
        } else {
          if (mu == 0) {
            if (check_inv(S_x(as.matrix(x_scaled[y_elems, ] / (n * h^d)), q, kernel_type))[1] == TRUE) {
              sx_mat <- solve(S_x(as.matrix(x_scaled[y_elems, ] / (n * h^d)), q, kernel_type))
            } else {
              singular_flag <- TRUE
              sx_mat <- matrix(0L, nrow = length(e_nu), ncol = length(e_nu))
            }
            bx <- b_x(as.matrix(x_scaled), sx_mat, e_nu, q, kernel_type)
            # sy matrix
            if (check_inv(S_x(as.matrix(y_scaled[y_elems]), p, kernel_type) / (n * h))[1] == TRUE) {
              sy_mat <- solve(S_x(as.matrix(y_scaled[y_elems]), p, kernel_type) / (n * h))
            } else {
              singular_flag <- TRUE
              sy_mat <- matrix(0L, nrow = length(e_mu), ncol = length(e_mu))
            }
            if (check_inv(S_x(as.matrix(yp_scaled[yp_elems]), p, kernel_type) / (n * h))[1] == TRUE) {
              syp_mat <- solve(S_x(as.matrix(yp_scaled[yp_elems]), p, kernel_type) / (n * h))
            } else {
              singular_flag <- TRUE
              syp_mat <- matrix(0L, nrow = length(e_mu), ncol = length(e_mu))
            }

            # computing y and yprime vectors
            a_y <- b_x(matrix(y_scaled[elems]), sy_mat, e_mu, p, kernel_type)
            a_yp <- b_x(matrix(yp_scaled[elems]), syp_mat, e_mu, p, kernel_type)

            # part 1 of the sum
            aj <- cumsum(a_y)
            ak <- cumsum(a_yp)
          } else {
            # sy matrix
            if (check_inv(S_x(as.matrix(y_scaled[y_elems]), p, kernel_type) / (n * h))[1] == TRUE) {
              sy_mat <- solve(S_x(as.matrix(y_scaled[y_elems]), p, kernel_type) / (n * h))
            } else {
              singular_flag <- TRUE
              sy_mat <- matrix(0L, nrow = length(e_mu), ncol = length(e_mu))
            }
            if (check_inv(S_x(as.matrix(yp_scaled[yp_elems]), p, kernel_type) / (n * h))[1] == TRUE) {
              syp_mat <- solve(S_x(as.matrix(yp_scaled[yp_elems]), p, kernel_type) / (n * h))
            } else {
              singular_flag <- TRUE
              syp_mat <- matrix(0L, nrow = length(e_mu), ncol = length(e_mu))
            }

            # computing y and yprime vectors
            a_y <- b_x(as.matrix(y_scaled[y_elems]), sy_mat, e_mu, p, kernel_type)
            a_yp <- b_x(as.matrix(yp_scaled[elems]), syp_mat, e_mu, p, kernel_type)

            # part 1 of the sum
            aj <- cumsum(a_y)
            ak <- cumsum(a_yp)
          }

          # populating matrix
          for (k in 1:length(elems)) {
            t_1 <- aj[k] * ak[k]
            t_2 <- (n - k) * a_yp[k] * aj[k]
            t_3 <- (n - k) * a_y[k] * ak[k]
            t_4 <- (n - k)^2 * a_y[k] * a_yp[k]
            c_hat[i, i] <- c_hat[i, i] + bx[1, elems[k]]^2 * (t_1 + t_2 + t_3 + t_4)
          }
        }
        # estimated means
        if (mu == 0) {
          theta_y <- fhat(
            x_data = x_data, y_data = y_data, x = x, y_grid = y, p = 2,
            q = 1, mu = 0, nu = 0, h = h, kernel_type = kernel_type
          )$est
          theta_yp <- fhat(
            x_data = x_data, y_data = y_data, x = x, y_grid = y_prime,
            p = 2, q = 1, mu = 0, nu = 0, h = h, kernel_type = kernel_type
          )$est
        } else {
          theta_y <- fhat(
            x_data = x_data, y_data = y_data, x = x, y_grid = y, p = p,
            q = q, mu = mu, nu = nu, h = h, kernel_type = kernel_type
          )$est
          theta_yp <- fhat(
            x_data = x_data, y_data = y_data, x = x, y_grid = y_prime,
            p = p, q = q, mu = mu, nu = nu, h = h, kernel_type = kernel_type
          )$est
        }

        # filling matrix, using symmetry
        c_hat[i, i] <- c_hat[i, i] / (n * (n - 1)^2) - theta_y * theta_yp / n^2
      }
    } else if (cov_flag == "full") {
      for (i in 1:ng) {
        for (j in 1:i) {
          # relevant entries w.r.t. y and y_prime
          y <- y_grid[i]
          y_prime <- y_grid[j]
          y_scaled <- (y_sorted - y) / h
          yp_scaled <- (y_sorted - y_prime) / h
          y_elems <- which(abs(y_scaled) <= 1)
          yp_elems <- which(abs(yp_scaled) <= 1)
          elems <- intersect(y_elems, yp_elems)

          if (length(elems) <= 5) {
            c_hat[i, j] <- 0
            c_hat[j, i] <- 0
          } else {
            if (mu == 0) {
              if (check_inv(S_x(as.matrix(x_scaled[y_elems, ] / (n * h^d)), q, kernel_type))[1] == TRUE) {
                sx_mat <- solve(S_x(as.matrix(x_scaled[y_elems, ] / (n * h^d)), q, kernel_type))
              } else {
                singular_flag <- TRUE
                sx_mat <- matrix(0L, nrow = length(e_nu), ncol = length(e_nu))
              }
              bx <- b_x(as.matrix(x_scaled), sx_mat, e_nu, q, kernel_type)
              # sy matrix
              if (check_inv(S_x(as.matrix(y_scaled[y_elems]), p, kernel_type) / (n * h))[1] == TRUE) {
                sy_mat <- solve(S_x(as.matrix(y_scaled[y_elems]), p, kernel_type) / (n * h))
              } else {
                singular_flag <- TRUE
                sy_mat <- matrix(0L, nrow = length(e_mu), ncol = length(e_mu))
              }
              if (check_inv(S_x(as.matrix(yp_scaled[yp_elems]), p, kernel_type) / (n * h))[1] == TRUE) {
                syp_mat <- solve(S_x(as.matrix(yp_scaled[yp_elems]), p, kernel_type) / (n * h))
              } else {
                singular_flag <- TRUE
                syp_mat <- matrix(0L, nrow = length(e_mu), ncol = length(e_mu))
              }

              # computing y and yprime vectors
              a_y <- b_x(matrix(y_scaled[elems]), sy_mat, e_mu, p, kernel_type)
              a_yp <- b_x(matrix(yp_scaled[elems]), syp_mat, e_mu, p, kernel_type)

              # part 1 of the sum
              aj <- cumsum(a_y)
              ak <- cumsum(a_yp)
            } else {
              # sy matrix
              if (check_inv(S_x(as.matrix(y_scaled[y_elems]), p, kernel_type) / (n * h))[1] == TRUE) {
                sy_mat <- solve(S_x(as.matrix(y_scaled[y_elems]), p, kernel_type) / (n * h))
              } else {
                singular_flag <- TRUE
                sy_mat <- matrix(0L, nrow = length(e_mu), ncol = length(e_mu))
              }
              if (check_inv(S_x(as.matrix(yp_scaled[yp_elems]), p, kernel_type) / (n * h))[1] == TRUE) {
                syp_mat <- solve(S_x(as.matrix(yp_scaled[yp_elems]), p, kernel_type) / (n * h))
              } else {
                singular_flag <- TRUE
                syp_mat <- matrix(0L, nrow = length(e_mu), ncol = length(e_mu))
              }

              # computing y and yprime vectors
              a_y <- b_x(as.matrix(y_scaled[y_elems]), sy_mat, e_mu, p, kernel_type)
              a_yp <- b_x(as.matrix(yp_scaled[elems]), syp_mat, e_mu, p, kernel_type)

              # part 1 of the sum
              aj <- cumsum(a_y)
              ak <- cumsum(a_yp)
            }

            # populating matrix
            for (k in 1:length(elems)) {
              t_1 <- aj[k] * ak[k]
              t_2 <- (n - k) * a_yp[k] * aj[k]
              t_3 <- (n - k) * a_y[k] * ak[k]
              t_4 <- (n - k)^2 * a_y[k] * a_yp[k]
              c_hat[i, j] <- c_hat[i, j] + bx[1, elems[k]]^2 * (t_1 + t_2 + t_3 + t_4)
            }
          }
          # estimated means
          if (mu == 0) {
            theta_y <- fhat(
              x_data = x_data, y_data = y_data, x = x, y_grid = y, p = 2,
              q = 1, mu = 0, nu = 0, h = h, kernel_type = kernel_type
            )$est
            theta_yp <- fhat(
              x_data = x_data, y_data = y_data, x = x, y_grid = y_prime,
              p = 2, q = 1, mu = 0, nu = 0, h = h, kernel_type = kernel_type
            )$est
          } else {
            theta_y <- fhat(
              x_data = x_data, y_data = y_data, x = x, y_grid = y, p = p,
              q = q, mu = mu, nu = nu, h = h, kernel_type = kernel_type
            )$est
            theta_yp <- fhat(
              x_data = x_data, y_data = y_data, x = x, y_grid = y_prime,
              p = p, q = q, mu = mu, nu = nu, h = h, kernel_type = kernel_type
            )$est
          }

          # filling matrix, using symmetry
          c_hat[i, j] <- c_hat[i, j] / (n * (n - 1)^2) - theta_y * theta_yp / n^2
          c_hat[j, i] <- c_hat[i, j]
        }
      }
    }
    if (mu == 0) {
      c_hat <- sweep(sweep(c_hat, MARGIN = 1, FUN = "*", STATS = 1 / (n * h^(d + mu + nu))),
        MARGIN = 2, FUN = "*", STATS = 1 / (n * h^(d + mu + nu))
      )
    } else {
      c_hat <- sweep(sweep(c_hat, MARGIN = 1, FUN = "*", STATS = 1 / (n * h^(d + mu + nu + 1))),
        MARGIN = 2, FUN = "*", STATS = 1 / (n * h^(d + mu + nu + 1))
      )
      diag(c_hat) <- diag(c_hat / 2)
    }
  } else {
    # initialize matrix
    c_hat <- matrix(0L, nrow = ng, ncol = ng)

    for (i in 1:ng) {
      for (j in 1:i) {
        hx <- mean(h[i], h[j])
        # localization for x
        idx <- which(rowSums(abs(sweep(x_data, 2, x)) <= hx) == d)

        x_data_loc <- matrix(x_data[idx, ], ncol = d)
        y_data_loc <- y_data[idx]

        # idx of ordering wrt y
        sort_idx <- sort(y_data_loc, index.return = TRUE)$ix

        # sorting datasets
        y_data_loc <- matrix(y_data_loc[sort_idx])
        x_data_loc <- matrix(x_data_loc[sort_idx, ], ncol = d)

        # x constants
        x_scaled <- sweep(x_data_loc, 2, x) / (hx^d)
        if (check_inv(S_x(matrix(x_scaled, ncol = d), q, kernel_type) / (n * hx^d))[1] == TRUE) {
          sx_mat <- solve(S_x(matrix(x_scaled, ncol = d), q, kernel_type) / (n * hx^d))
        } else {
          singular_flag <- TRUE
          sx_mat <- matrix(0L, nrow = length(e_nu), ncol = length(e_nu))
        }
        bx <- b_x(matrix(x_scaled, ncol = d), sx_mat, e_nu, q, kernel_type)

        # relevant entries wrt y and y_prime
        y <- y_grid[i]
        y_prime <- y_grid[j]
        y_scaled <- (y_data_loc - y) / h[i]
        yp_scaled <- (y_data_loc - y_prime) / h[j]
        y_elems <- which(abs(y_scaled) <= 1)
        yp_elems <- which(abs(yp_scaled) <= 1)
        elems <- intersect(y_elems, yp_elems)

        if (length(elems) <= 5) {
          c_hat[i, j] <- 0
          c_hat[j, i] <- 0
        } else {
          if (mu == 0) {
            if (check_inv(S_x(as.matrix(x_scaled[y_elems, ] / (n * hx^d)), q, kernel_type))[1] == TRUE) {
              sx_mat <- solve(S_x(as.matrix(x_scaled[y_elems, ] / (n * hx^d)), q, kernel_type))
            } else {
              singular_flag <- TRUE
              sx_mat <- matrix(0L, nrow = length(e_nu), ncol = length(e_nu))
            }
            bx <- b_x(as.matrix(x_scaled), sx_mat, e_nu, q, kernel_type)
            # sy matrix
            if (check_inv(S_x(as.matrix(y_scaled[y_elems]), p, kernel_type) / (n * h[i]))[1] == TRUE) {
              sy_mat <- solve(S_x(as.matrix(y_scaled[y_elems]), p, kernel_type) / (n * h[i]))
            } else {
              singular_flag <- TRUE
              sy_mat <- matrix(0L, nrow = length(e_mu), ncol = length(e_mu))
            }
            if (check_inv(S_x(as.matrix(yp_scaled[yp_elems]), p, kernel_type) / (n * h[j]))[1] == TRUE) {
              syp_mat <- solve(S_x(as.matrix(yp_scaled[yp_elems]), p, kernel_type) / (n * h[j]))
            } else {
              singular_flag <- TRUE
              syp_mat <- matrix(0L, nrow = length(e_mu), ncol = length(e_mu))
            }

            # computing y and yprime vectors
            a_y <- b_x(matrix(y_scaled[elems]), sy_mat, e_mu, p, kernel_type)
            a_yp <- b_x(matrix(yp_scaled[elems]), syp_mat, e_mu, p, kernel_type)

            # part 1 of the sum
            aj <- cumsum(a_y)
            ak <- cumsum(a_yp)
          } else {
            # sy matrix
            if (check_inv(S_x(as.matrix(y_scaled[y_elems]), p, kernel_type) / (n * h[i]))[1] == TRUE) {
              sy_mat <- solve(S_x(as.matrix(y_scaled[y_elems]), p, kernel_type) / (n * h[i]))
            } else {
              singular_flag <- TRUE
              sy_mat <- matrix(0L, nrow = length(e_mu), ncol = length(e_mu))
            }
            if (check_inv(S_x(as.matrix(yp_scaled[yp_elems]), p, kernel_type) / (n * h[j]))[1] == TRUE) {
              syp_mat <- solve(S_x(as.matrix(yp_scaled[yp_elems]), p, kernel_type) / (n * h[j]))
            } else {
              singular_flag <- TRUE
              syp_mat <- matrix(0L, nrow = length(e_mu), ncol = length(e_mu))
            }

            # computing y and yprime vectors
            a_y <- b_x(matrix(y_scaled[elems]), sy_mat, e_mu, p, kernel_type)
            a_yp <- b_x(matrix(yp_scaled[elems]), syp_mat, e_mu, p, kernel_type)

            # part 1 of the sum
            aj <- cumsum(a_y)
            ak <- cumsum(a_yp)
          }

          # populating matrix
          for (k in 1:length(elems)) {
            t_1 <- aj[k] * ak[k]
            t_2 <- (n - k) * a_yp[k] * aj[k]
            t_3 <- (n - k) * a_y[k] * ak[k]
            t_4 <- (n - k)^2 * a_y[k] * a_yp[k]
            c_hat[i, j] <- c_hat[i, j] + bx[1, elems[k]]^2 * (t_1 + t_2 + t_3 + t_4)
          }
        }
        # estimated means
        if (mu == 0) {
          theta_y <- fhat(
            x_data = x_data, y_data = y_data, x = x, y_grid = y, p = 2,
            q = 1, mu = 0, nu = 0, h = h[i], kernel_type = kernel_type
          )$est
          theta_yp <- fhat(
            x_data = x_data, y_data = y_data, x = x, y_grid = y_prime,
            p = 2, q = 1, mu = 0, nu = 0, h = h[j], kernel_type = kernel_type
          )$est
        } else {
          theta_y <- fhat(
            x_data = x_data, y_data = y_data, x = x, y_grid = y, p = p,
            q = q, mu = mu, nu = nu, h = h[i], kernel_type = kernel_type
          )$est
          theta_yp <- fhat(
            x_data = x_data, y_data = y_data, x = x, y_grid = y_prime,
            p = p, q = q, mu = mu, nu = nu, h = h[j], kernel_type = kernel_type
          )$est
        }

        # filling matrix, using symmetry
        c_hat[i, j] <- c_hat[i, j] / (n * (n - 1)^2) - theta_y * theta_yp / n^2
        c_hat[j, i] <- c_hat[i, j]
      }
    }

    if (mu == 0) {
      c_hat <- sweep(sweep(c_hat, MARGIN = 1, FUN = "*", STATS = 1 / (n * h^(d + mu + nu))),
        MARGIN = 2, FUN = "*", STATS = 1 / (n * h^(d + mu + nu))
      )
    } else {
      c_hat <- sweep(sweep(c_hat, MARGIN = 1, FUN = "*", STATS = 1 / (n * h^(d + mu + nu + 1))),
        MARGIN = 2, FUN = "*", STATS = 1 / (n * h^(d + mu + nu + 1))
      )
      diag(c_hat) <- diag(c_hat / 2)
    }
  }

  return(list("cov" = c_hat, "singular_flag" = singular_flag))
}
#######################################################################################
# Supplemental Functions
#######################################################################################

#' @title bx
#' @description Function for estimating the constants in the estimation formula
#' @param datavec Dataset, vector.
#' @param s_mat S_hat matrix.
#' @param q Polynomial order.
#' @param kernel_type Kernel function choice.
#' @return Vector of products for each data point.
#' @keywords internal
b_x <- function(datavec, s_mat, e_vec, q, kernel_type) {
  eff_n <- length(datavec[, 1])
  Rq <- matrix(0L, ncol = eff_n)
  for (i in 1:eff_n) {
    Rq[i] <- (poly_base(datavec[i, ], q) * kernel_eval(datavec[i, ], kernel_type)) %*% (s_mat %*% e_vec)
  }
  return(Rq)
}
