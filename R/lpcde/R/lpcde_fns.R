#######################################################################################
# This file contains code for generating point estimates (Internal Functions)
#######################################################################################

#' @title lpcde_fn: conditional density estimator.
#' @description Function for estimating the density function and its derivatives.
#' @param y_data response variable dataset, vector.
#' @param x_data covariate dataset, vector or matrix.
#' @param y_grid Numeric vector, specifies the grid of evaluation points along y-direction.
#' @param x Numeric vector or matrix, specifies the grid of evaluation points along x-direction.
#' @param p polynomial order for y.
#' @param q polynomial order for covariates.
#' @param p_RBC Nonnegative integer, specifies the order of the local polynomial for \code{Y} used to
#' construct bias-corrected point estimates. (Default is \code{p+1}.)
#' @param q_RBC Nonnegative integer, specifies the order of the local polynomial for \code{X} used to
#' construct bias-corrected point estimates. (Default is \code{q+1}.)
#' @param bw Numeric, bandwidth vector.
#' @param mu degree of derivative with respect to y.
#' @param nu degree of derivative with respect to x.
#' @param kernel_type kernel function choice.
#' @param rbc Boolean for whether to return RBC estimate and standard errors.
# @param var_type String. type of variance estimator to implement.
# choose from "ustat" and "asymp".
#' @return conditional density estimate at all grid points.
#' @keywords internal
lpcde_fn = function(y_data, x_data, y_grid, x, p, q, p_RBC, q_RBC, bw, mu, nu,
                    kernel_type, rbc = FALSE){
  # initializing output vectors
  est = matrix(0L, nrow = length(y_grid), ncol = 1)
  se = matrix(0L, nrow = length(y_grid), ncol = 1)

  # estimate
  f_hat_val = fhat(x_data=x_data, y_data=y_data, x=x, y_grid=y_grid, p=p, q=q,
                   mu=mu, nu=nu, h=bw, kernel_type=kernel_type)
  est = f_hat_val$est
  eff.n = f_hat_val$eff.n

  # covariance matrix

  # standard errors
  # if(var_type=="ustat"){
    covMat = cov_hat(x_data=x_data, y_data=y_data, x=x, y_grid=y_grid, p=p, q=q,
                   mu=mu, nu=nu, h=bw, kernel_type=kernel_type)
    se = sqrt(abs(diag(covMat)))
  # }else if (var_type == "simul"){
  #   covMat = simul_cov(x_data=x_data, y_data=y_data, x=x, y_grid=y_grid, p=p, q=q,
  #                  mu=mu, nu=nu, h=bw, kernel_type=kernel_type)
  #   se = sqrt(abs(diag(covMat$covMat)))
  # }else if (var_type == "asymp"){
    # covMat = asymp_var(x_data=x_data, y_data=y_data, x=x, y_grid=y_grid, p=p, q=q,
                       # mu=mu, nu=nu, h=bw, kernel_type=kernel_type)
    # se = sqrt(abs(diag(covMat)))
  # }

  if (rbc){
    est_rbc = matrix(0L, nrow = length(y_grid), ncol = 1)
    se_rbc = matrix(0L, nrow = length(y_grid), ncol = 1)
    if(p_RBC == p && q_RBC==q){
      est_rbc = est
      se_rbc = se
      covMat_rbc = covMat
    }else{
      # estimate
      est_rbc = fhat(x_data=x_data, y_data=y_data, x=x, y_grid=y_grid, p=p_RBC,
                     q=q_RBC, mu=mu, nu=nu, h=bw, kernel_type=kernel_type)$est

      # covariance matrix

      # standard errors
      # if(var_type=="ustat"){
        covMat_rbc = cov_hat(x_data=x_data, y_data=y_data, x=x, y_grid=y_grid,
                           p=p_RBC, q=q_RBC, mu=mu, nu=nu, h=bw,
                           kernel_type=kernel_type)
        se_rbc = sqrt(abs(diag(covMat_rbc)))
      # }else if (var_type == "simul"){
      #   covMat_rbc = simul_cov(x_data=x_data, y_data=y_data, x=x, y_grid=y_grid, p=p_RBC, q=q_RBC,
      #                  mu=mu, nu=nu, h=bw, kernel_type=kernel_type)
      #   se_rbc = sqrt(abs(diag(covMat$covMat)))
      # }else if (var_type == "asymp"){
        # covMat_rbc = asymp_var(x_data=x_data, y_data=y_data, x=x, y_grid=y_grid, p=p_RBC, q=q_RBC,
                           # mu=mu, nu=nu, h=bw, kernel_type=kernel_type)
        # se_rbc = sqrt(abs(diag(covMat_rbc)))
      # }
    }

    # with rbc results
    estimate = cbind(y_grid, bw, est, est_rbc, se, se_rbc)
    colnames(estimate) = c("y_grid","bw", "est","est_RBC", "se", "se_RBC")
    rownames(estimate) = c()
    est_result = list("est" = estimate,
                      "CovMat" = list("CovMat" = covMat,
                                      "CovMat_RBC" = covMat_rbc),
                      "x" = x, "eff_n" = eff.n)
  } else{
    # without rbc results
    # setting rbc values to be the same as non-rbc values
    est_rbc = est
    se_rbc = se
    covMat_rbc = covMat
    # generating matrices and list to return
    estimate = cbind(y_grid, bw, est, est_rbc, se, se_rbc)
    colnames(estimate) = c("y_grid","bw", "est","est_RBC", "se", "se_RBC")
    rownames(estimate) = c()
    est_result = list("est" = estimate,
                      "CovMat" = list("CovMat" = covMat,
                                      "CovMat_RBC" = covMat_rbc),
                      "x" = x, "eff_n" = eff.n)
  }

  return(est_result)
}

#' @title fhat: estimator
#' @description Function for estimating the density function and its derivatives.
#' @param y_data response variable dataset, vector.
#' @param x_data covariate dataset, vector or matrix.
#' @param y_grid Numeric vector, specifies the grid of evaluation points along y-direction.
#' @param x Numeric vector or matrix, specifies the grid of evaluation points along x-direction.
#' @param p polynomial order for y.
#' @param q polynomial order for covariates.
#' @param h Numeric, bandwidth vector.
#' @param mu degree of derivative with respect to y.
#' @param nu degree of derivative with respect to x.
#' @param kernel_type kernel function choice.
#' @return conditional density estimate at all grid points.
#' @keywords internal
fhat = function(x_data, y_data, x, y_grid, p, q, mu, nu, h, kernel_type){
  # setting constants
  n = length(y_data)
  d = ncol(x_data)
  ng = length(y_grid)

  # x basis vector
  e_nu = basis_vec(x, q, nu)

  # y basis vector
  e_mu = basis_vec(0, p, mu)

  f_hat = matrix(0L, nrow = ng)
  nh_vec = matrix(0L, nrow = ng)

  if (length(unique(h))==1){
    h = h[1]
    # localization for x
    idx = which(abs(x_data-x)<=h^d)

    x_idx = x_data[idx]
    y_idx = y_data[idx]

    # idx of ordering wrt y
    sort_idx = sort(y_idx, index.return=TRUE)$ix

    # sorting datasets
    y_sorted = as.matrix(y_idx[sort_idx])
    x_sorted = as.matrix(x_idx[sort_idx])

    # x constants
    x_scaled = (x_sorted-x)/(h^d)
    sx_mat = solve(S_x(x_scaled, q, kernel_type)/(n*h^d))
    bx = b_x(x_scaled, sx_mat, e_nu, q, kernel_type)

    for (j in 1:ng){
      y = y_grid[j]
      y_scaled = (y_sorted-y)/h
      y_elems = which(abs(y_scaled)<=1)

      if (length(y_elems)<=5){
        f_hat[j] = 0
        nh_vec [j] = length(y_elems)
      } else {
        if (mu==0){
          sx_mat = solve(S_x(as.matrix(x_scaled[y_elems]), q, kernel_type)/(n*h^d))
          bx = b_x(x_scaled, sx_mat, e_nu, q, kernel_type)
          # sy matrix
          sy_mat = solve(S_x(as.matrix(y_scaled), p, kernel_type)/(n*h))

          # y constants
          ax = b_x(y_scaled, sy_mat, e_mu, p, kernel_type)

          # adding and multiplying
          f_hat[j] = ax[y_elems] %*% cumsum(bx[y_elems])

          # number of datapoints used
          nh_vec[j] = length(y_elems)
        } else {
          # sy matrix
          sy_mat = solve(S_x(as.matrix(y_scaled[y_elems]), p, kernel_type)/(n*h))

          # y constants
          ax = b_x(y_scaled[y_elems], sy_mat, e_mu, p, kernel_type)

          # adding and multiplying
          f_hat[j] = ax %*% cumsum(bx[y_elems])

          # number of datapoints used
          nh_vec[j] = length(y_elems)
        }
      }
    }

    # scaling
    f_hat = f_hat/(n^2*h^(d+mu+nu+1))
  } else {
    for (j in 1:ng){
      # localization for x
      idx = which(abs(x_data-x)<=h[j])

      x_data_loc = x_data[idx]
      y_data_loc = y_data[idx]

      # idx of ordering wrt y
      sort_idx = sort(y_data_loc, index.return=TRUE)$ix

      # sorting datasets
      y_data_loc = as.matrix(y_data_loc[sort_idx])
      x_data_loc = as.matrix(x_data_loc[sort_idx])

      # x constants
      x_scaled = (x_data_loc-x)/(h[j]^d)
      sx_mat = solve(S_x(as.matrix(x_scaled), q, kernel_type)/(n*h[j]^d))
      bx = b_x(as.matrix(x_scaled), sx_mat, e_nu, q, kernel_type)

      y = y_grid[j]
      y_scaled = (y_data_loc-y)/h[j]
      y_elems = which(abs(y_scaled)<=1)

      # sy matrix
      sy_mat = solve(S_x(as.matrix(y_scaled[y_elems]), p, kernel_type)/(n*h[j]))

      # y constants
      ax = b_x(y_scaled[y_elems], sy_mat, e_mu, p, kernel_type)

      # adding and multiplying
      f_hat[j] = ax %*% cumsum(bx[y_elems])

      # scaling
      f_hat[j] = f_hat[j]/(n^2*h[j]^(d+mu+nu+1))
      # number of datapoints used
      nh_vec[j] = length(y_elems)
    }
  }
  return(list("est" = f_hat, "eff.n" = nh_vec))

}

##' @title asymp_var: asymptotic covariance estimator
##' @description Function for estimating the variance-covariance matrix.
##' @param y_data response variable dataset, vector.
##' @param x_data covariate dataset, vector or matrix.
##' @param y_grid Numeric vector, specifies the grid of evaluation points along y-direction.
##' @param x Numeric vector or matrix, specifies the grid of evaluation points along x-direction.
##' @param p polynomial order for y.
##' @param q polynomial order for covariates.
##' @param h Numeric, bandwidth vector.
##' @param mu degree of derivative with respect to y.
##' @param nu degree of derivative with respect to x.
##' @param kernel_type kernel function choice.
##' @return covariance matrix for all the grid points
##' @keywords internal
#asymp_var = function(x_data, y_data, x, y_grid, p, q, mu, nu, h, kernel_type){

#  n = length(y_data)
#  d = ncol(x_data)
#  ng = length(y_grid)

#  # x basis vector
#  e_nu = basis_vec(x, q, nu)

#  # y basis vector
#  e_mu = basis_vec(1, p, mu)

#  if(length(unique(h))==1){
#    h = h[1]
#    theta = fhat(x_data=as.matrix(x_data), y_data=as.matrix(y_data), x=x,
#                 y_grid=y_grid, p=2, q=1, mu=1, nu=0, h=h,
#                 kernel_type=kernel_type)$est
#    theta_00 = fhat(x_data=as.matrix(x_data), y_data=as.matrix(y_data), x=x,
#                    y_grid=y_grid, p=1, q=1, mu=0, nu=0, h=h,
#                    kernel_type=kernel_type)$est

#    # localization for x
#    idx = which(abs(x_data-x)<=h)

#    x_data = x_data[idx]
#    y_data = y_data[idx]

#    # idx of ordering wrt y
#    sort_idx = sort(y_data, index.return=TRUE)$ix

#    # sorting datasets
#    y_data = as.matrix(y_data[sort_idx])
#    x_data = as.matrix(x_data[sort_idx])

#    # x constants
#    x_scaled = (x_data-x)/(h^d)

#    sx_mat = solve(S_x(as.matrix(x_scaled), q, kernel_type)/(n*h^d))
#    Tx = T_x(x_scaled, x, q, h, kernel_type)/h^d

#    c_mat = matrix(0L, nrow=ng, ncol=ng)
#    for (j in 1:ng){
#      for (k in 1:j){
#        y = y_grid[j]
#        yp = y_grid[k]
#        y_scaled = (y_data-y)/h
#        yp_scaled = (y_data-yp)/h
#        y_elems = which(abs(y_scaled)<1)
#        yp_elems = which(abs(yp_scaled)<1)
#        sy_mat = solve(S_x(as.matrix(y_scaled), p, kernel_type)/(n*h))
#        syp_mat = solve(S_x(as.matrix(yp_scaled), p, kernel_type)/(n*h))
#        elems = intersect(y_elems, yp_elems)
#        if (length(elems)==0){
#          c_mat[j, k] = 0
#          c_mat[k, j] = c_mat[j, k]
#        }else{
#          if (mu==0){
#            c_mat[j, k] = theta_00[j]*(1-theta_00[j]) * (sx_mat%*%Tx%*%sx_mat)[nu+1, nu+1]
#            c_mat[k, j] = c_mat[j, k]
#          }else{
#            Ty = T_y(y_scaled[elems], yp_scaled[elems], p, kernel_type)/(choose(n, 2)*h^2)
#            c_mat[j, k] = theta[j] * (sy_mat%*%Ty%*%syp_mat)[mu+1, mu+1] * (sx_mat%*%Tx%*%sx_mat)[nu+1, nu+1]
#            c_mat[k, j] = c_mat[j, k]
#          }
#        }
#      }
#    }
#    if(mu == 0){
#      c_mat = sweep(sweep(c_mat, MARGIN=1, FUN="*", STATS=1/(n*h^(d+2*nu))),
#                                            MARGIN=2, FUN="*", STATS=1/(n*h^(d+2*nu)))
#    }else{
#      c_mat = sweep(sweep(c_mat, MARGIN=1, FUN="*", STATS=1/(n*h^(d+2*mu+2*nu-1))),
#                    MARGIN=2, FUN="*", STATS=1/(n*h^(d+2*mu+2*nu-1)))
#    }
#  }else{

#    theta = fhat(x_data=as.matrix(x_data), y_data=as.matrix(y_data), x=x, y_grid=y_grid, p=2,
#                      q=1, mu=1, nu=0, h=h, kernel_type=kernel_type)$est
#    theta_00 = fhat(x_data=as.matrix(x_data), y_data=as.matrix(y_data), x=x, y_grid=y_grid, p=3,
#                            q=1, mu=0, nu=0, h=h, kernel_type=kernel_type)$est
#    # localization for x
#    idx = which(abs(x_data-x)<=mean(h))

#    x_data = x_data[idx]
#    y_data = y_data[idx]

#    # idx of ordering wrt y
#    sort_idx = sort(y_data, index.return=TRUE)$ix

#    # sorting datasets
#    y_data = as.matrix(y_data[sort_idx])
#    x_data = as.matrix(x_data[sort_idx])

#    # x constants
#    x_scaled = (x_data-x)/(mean(h)^d)

#    sx_mat = solve(S_x(as.matrix(x_scaled), q, kernel_type)/(n*mean(h)^d))
#    Tx = T_x(x_scaled, x, q, mean(h), kernel_type)/(mean(h)^d)

#    c_mat = matrix(0L, nrow=ng, ncol=ng)
#    for (j in 1:ng){
#      for (k in 1:j){
#        y = y_grid[j]
#        yp = y_grid[k]
#        y_scaled = (y_data-y)/h[j]
#        yp_scaled = (y_data-yp)/h[k]
#        y_elems = which(abs(y_scaled)<1)
#        yp_elems = which(abs(yp_scaled)<1)
#        sy_mat = solve(S_x(as.matrix(y_scaled), p, kernel_type)/(n*h[j]))
#        syp_mat = solve(S_x(as.matrix(yp_scaled), p, kernel_type)/(n*h[k]))
#        elems = intersect(y_elems, yp_elems)
#        if (length(elems)==0){
#          c_mat[j, k] = 0
#          c_mat[k, j] = c_mat[j, k]
#        }else {
#          if (mu==0){
#            c_mat[j, k] = theta_00[j]*(1-theta_00[j]) * (sx_mat%*%Tx%*%sx_mat)[nu+1, nu+1]
#            c_mat[k, j] = c_mat[j, k]
#          }else{
#            Ty = T_y(y_scaled[elems], yp_scaled[elems], p, kernel_type)/(choose(n, 2)*h[j]*h[k])
#            c_mat[j, k] = theta[j] * (sy_mat%*%Ty%*%syp_mat)[mu+1, mu+1] * (sx_mat%*%Tx%*%sx_mat)[nu+1, nu+1]
#            c_mat[k, j] = c_mat[j, k]
#          }
#        }
#      }
#    }
#    if(mu == 0){
#      c_mat = sweep(sweep(c_mat, MARGIN=1, FUN="*", STATS=1/(n*h^(d+2*nu))),
#                                            MARGIN=2, FUN="*", STATS=1/(n*h^(d+2*nu)))
#    }else{
#      c_mat = sweep(sweep(c_mat, MARGIN=1, FUN="*", STATS=1/(n*h^(d+2*mu+2*nu-1))),
#                    MARGIN=2, FUN="*", STATS=1/(n*h^(d+2*mu+2*nu-1)))
#    }
#  }
#  return(c_mat)
#}

#' @title cov_hat: covariance estimator
#' @description Function for estimating the variance-covariance matrix.
#' @param y_data response variable dataset, vector.
#' @param x_data covariate dataset, vector or matrix.
#' @param y_grid Numeric vector, specifies the grid of evaluation points along
#' y-direction.
#' @param x Numeric vector or matrix, specifies the grid of evaluation points
#' along x-direction.
#' @param p polynomial order for y.
#' @param q polynomial order for covariates.
#' @param h Numeric, bandwidth vector.
#' @param mu degree of derivative with respect to y.
#' @param nu degree of derivative with respect to x.
#' @param kernel_type kernel function choice.
#' @return covariance matrix for all the grid points
#' @keywords internal
cov_hat = function(x_data, y_data, x, y_grid, p, q, mu, nu, h, kernel_type){
  # setting constants
  n = length(y_data)
  d = ncol(x_data)
  ng = length(y_grid)

  # x basis vector
  e_nu = basis_vec(x, q, nu)

  # y basis vector
  e_mu = basis_vec(1, p, mu)

  if(length(unique(h))==1){
    h = h[1]
    # localization for x
    idx = which(abs(x_data-x)<=h)

    x_data = x_data[idx]
    y_data = y_data[idx]

    # idx of ordering wrt y
    sort_idx = sort(y_data, index.return=TRUE)$ix

    # sorting datasets
    y_data = as.matrix(y_data[sort_idx])
    x_data = as.matrix(x_data[sort_idx])

    # x constants
    x_scaled = (x_data-x)/(h^d)
    sx_mat = solve(S_x(as.matrix(x_scaled), q, kernel_type)/(n*h^d))
    bx = b_x(as.matrix(x_scaled), sx_mat, e_nu, q, kernel_type)

    # initialize matrix
    c_hat = matrix(0L, nrow=ng, ncol=ng)

    for (i in 1:ng){
      for (j in 1:i){
        # relevant entries wrt y and y_prime
        y = y_grid[i]
        y_prime = y_grid[j]
        y_scaled = (y_data - y)/h
        yp_scaled = (y_data-y_prime)/h
        y_elems = which(abs(y_scaled)<=1)
        yp_elems = which(abs(yp_scaled)<=1)
        elems = intersect(y_elems, yp_elems)

        if ( length(elems) <= 5){
          c_hat[i, j] = c_hat[i, j]
        } else{
          if (mu==0){
            sx_mat = solve(S_x(as.matrix(x_scaled[y_elems]), q, kernel_type)/(n*h^d))
            bx = b_x(as.matrix(x_scaled), sx_mat, e_nu, q, kernel_type)
            # sy matrix
            sy_mat = solve(S_x(as.matrix(y_scaled[y_elems]), p, kernel_type)/(n*h))
            syp_mat = solve(S_x(as.matrix(yp_scaled[yp_elems]), p, kernel_type)/(n*h))

            # computing y and yprime vectors
            a_y = b_x(y_scaled[elems], sy_mat, e_mu, p, kernel_type)
            a_yp =  b_x(yp_scaled[elems], syp_mat, e_mu, p, kernel_type)

            # part 1 of the sum
            aj = cumsum(a_y)
            ak = cumsum(a_yp)

          }else{
            # sy matrix
            sy_mat = solve(S_x(as.matrix(y_scaled[y_elems]), p, kernel_type)/(n*h))
            syp_mat = solve(S_x(as.matrix(yp_scaled[yp_elems]), p, kernel_type)/(n*h))

            # computing y and yprime vectors
            a_y = b_x(y_scaled[elems], sy_mat, e_mu, p, kernel_type)
            a_yp =  b_x(yp_scaled[elems], syp_mat, e_mu, p, kernel_type)

            # part 1 of the sum
            aj = cumsum(a_y)
            ak = cumsum(a_yp)
          }

          # populating matrix
          for (k in 1:length(elems)){
            t_1 = aj[k] * ak[k]
            t_2 = (n-k) * a_yp[k] * aj[k]
            t_3 = (n-k) * a_y[k] * ak[k]
            t_4 = (n-k)^2 * a_y[k] * a_yp[k]
            c_hat[i, j] = c_hat[i, j] + bx[elems[k]]^2 * (t_1 + t_2 + t_3 + t_4)
          }
        }
        # estimated means
        if (mu==0){
          theta_y = fhat(x_data=x_data, y_data=y_data, x=x, y_grid=y, p=2,
                         q=1, mu=0, nu=0, h=h, kernel_type=kernel_type)$est
          theta_yp = fhat(x_data=x_data, y_data=y_data, x=x, y_grid=y_prime,
                          p=2, q=1, mu=0, nu=0, h=h, kernel_type=kernel_type)$est

        }else{
          theta_y = fhat(x_data=x_data, y_data=y_data, x=x, y_grid=y, p=p,
                         q=q, mu=mu, nu=nu, h=h, kernel_type=kernel_type)$est
          theta_yp = fhat(x_data=x_data, y_data=y_data, x=x, y_grid=y_prime,
                          p=p, q=q, mu=mu, nu=nu, h=h, kernel_type=kernel_type)$est
        }

        # filling matrix, using symmetry
        c_hat[i, j] = c_hat[i, j]/(n*(n-1)^2) - theta_y*theta_yp/n^2
        c_hat[j, i] = c_hat[i, j]


      }
    }
    if(mu==0){
      c_hat = sweep(sweep(c_hat, MARGIN=1, FUN="*", STATS=1/(n*h^(d+mu+nu))),
                  MARGIN=2, FUN="*", STATS=1/(n*h^(d+mu+nu)))
    }else{
      c_hat = sweep(sweep(c_hat, MARGIN=1, FUN="*", STATS=1/(n*h^(d+mu+nu+1))),
                  MARGIN=2, FUN="*", STATS=1/(n*h^(d+mu+nu+1)))
      diag(c_hat) = diag(c_hat/2)
    }

  }else {
    # initialize matrix
    c_hat = matrix(0L, nrow=ng, ncol=ng)

    for (i in 1:ng){
      for (j in 1:i){
        # localization for x
        idx = which(abs(x_data-x)<=h[i])

        x_data_loc = x_data[idx]
        y_data_loc = y_data[idx]

        # idx of ordering wrt y
        sort_idx = sort(y_data_loc, index.return=TRUE)$ix

        # sorting datasets
        y_data_loc = as.matrix(y_data_loc[sort_idx])
        x_data_loc = as.matrix(x_data_loc[sort_idx])

        # x constants
        x_scaled = (x_data_loc-x)/(h[i]^d)
        sx_mat = solve(S_x(x_scaled, q, kernel_type)/(n*h[i]^d))
        bx = b_x(x_scaled, sx_mat, e_nu, q, kernel_type)

        # relevant entries wrt y and y_prime
        y = y_grid[i]
        y_prime = y_grid[j]
        y_scaled = (y_data_loc - y)/h[i]
        yp_scaled = (y_data_loc-y_prime)/h[j]
        y_elems = which(abs(y_scaled)<=1)
        yp_elems = which(abs(yp_scaled)<=1)
        elems = intersect(y_elems, yp_elems)

        if ( length(elems) == 0){
          c_hat[i, j] = c_hat[i, j]
        } else{
          if (mu==0){
            sx_mat = solve(S_x(as.matrix(x_scaled[y_elems]), q, kernel_type)/(n*h[i]^d))
            bx = b_x(as.matrix(x_scaled), sx_mat, e_nu, q, kernel_type)
            # sy matrix
            sy_mat = solve(S_x(as.matrix(y_scaled[y_elems]), p, kernel_type)/(n*h[i]))
            syp_mat = solve(S_x(as.matrix(yp_scaled[yp_elems]), p, kernel_type)/(n*h[j]))

            # computing y and yprime vectors
            a_y = b_x(y_scaled[elems], sy_mat, e_mu, p, kernel_type)
            a_yp =  b_x(yp_scaled[elems], syp_mat, e_mu, p, kernel_type)

            # part 1 of the sum
            aj = cumsum(a_y)
            ak = cumsum(a_yp)

          }else{
            # sy matrix
            sy_mat = solve(S_x(as.matrix(y_scaled[y_elems]), p, kernel_type)/(n*h[i]))
            syp_mat = solve(S_x(as.matrix(yp_scaled[yp_elems]), p, kernel_type)/(n*h[j]))

            # computing y and yprime vectors
            a_y = b_x(y_scaled[elems], sy_mat, e_mu, p, kernel_type)
            a_yp =  b_x(yp_scaled[elems], syp_mat, e_mu, p, kernel_type)

            # part 1 of the sum
            aj = cumsum(a_y)
            ak = cumsum(a_yp)
          }

          # populating matrix
          for (k in 1:length(elems)){
            t_1 = aj[k] * ak[k]
            t_2 = (n-k) * a_yp[k] * aj[k]
            t_3 = (n-k) * a_y[k] * ak[k]
            t_4 = (n-k)^2 * a_y[k] * a_yp[k]
            c_hat[i, j] = c_hat[i, j] + bx[elems[k]]^2 * (t_1 + t_2 + t_3 + t_4)
          }
        }
        # estimated means
        theta_y = fhat(x_data=x_data, y_data=y_data, x=x, y_grid=y, p=p, q=q,
                       mu=mu, nu=nu, h=h[i], kernel_type=kernel_type)$est
        theta_yp = fhat(x_data=x_data, y_data=y_data, x=x, y_grid=y_prime, p=p,
                        q=q, mu=mu, nu=nu, h=h[j], kernel_type=kernel_type)$est

        # filling matrix, using symmetry
        c_hat[i, j] = c_hat[i, j]/(n*(n-1)^2) - theta_y*theta_yp/n
        c_hat[j, i] = c_hat[i, j]
      }
    }
    # c_hat = c_hat/(n^2*h^(2*d+2*mu+2*nu+2))

    if(mu==0){
      c_hat = sweep(sweep(c_hat, MARGIN=1, FUN="*", STATS=1/(n*h^(d+mu+nu))),
                  MARGIN=2, FUN="*", STATS=1/(n*h^(d+mu+nu)))
    }else{
      c_hat = sweep(sweep(c_hat, MARGIN=1, FUN="*", STATS=1/(n*h^(d+mu+nu+1))),
                  MARGIN=2, FUN="*", STATS=1/(n*h^(d+mu+nu+1)))
      diag(c_hat) = diag(c_hat/2)
    }
  }

  return(c_hat)
}
#######################################################################################
# Supplemental Functions
#######################################################################################

#' @title bx = constants for each data and evaluation point pair
#' @description Function for estimating the constants in the estimation formula
#' @param data_vec dataset, vector.
#' @param s_mat s_hat matrix.
#' @param q polynomial order.
#' @param kernel_type kernel function choice.
#' @return vector of products for each data point.
#' @keywords internal
b_x = function(datavec, s_mat, e_vec, q, kernel_type){
  eff_n = length(datavec)
  Rq = matrix(0L, ncol = eff_n)
  for (i in 1:eff_n){
    Rq[i] = (poly_base(datavec[i], q) * kernel_eval(datavec[i], kernel_type)) %*% (s_mat %*% e_vec)
  }
  return(Rq)
}

# bias_corr = function(y_data, x_data, p, q, mu, nu, h, kernel_type){
#   y_data= as.matrix(y_data)
#   x_data = as.matrix(x_data)
#   n = length(y_data)
#   d=1
#   # x basis vector
#   e_nu = basis_vec(x, q, nu)
#
#   # y basis vector
#   e_mu = basis_vec(1, p, mu)
#   sy_mat = solve(S_x(y_data, p, kernel_type)/(n*h))
#   sx_mat = solve(S_x(x_data, q, kernel_type)/(n*h^d))
#   u = matrix(0L, nrow=n)
#   for(i in 1:n){
#     y_i = y_data[i]
#     x_i = x_data[i]
#     for (j in 1:n){
#       y_j = y_data[j]
#       if(j!=i){
#         u[i] = uij(y_i, y_j, x_i, y, x, kernel_type, symat, sxmat, p, q, mu, nu)
#       }
#     }
#   }
# }
