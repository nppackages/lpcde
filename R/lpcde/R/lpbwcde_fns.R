#######################################################################################
#' ROT Bandwidth selection
#'
#' Internal Function
#' @param y_data Numeric matrix/data frame, the raw data of independent.
#' @param x_data Numeric matrix/data frame, the raw data of covariates.
#' @param y_grid Numeric vector, the evaluation points.
#' @param x Numeric, specifies the evaluation point(s) in the x-direction.
#' @param p Integer, polynomial order.
#' @param q Integer, polynomial order.
#' @param mu Integer, order of derivative.
#' @param nu Integer, order of derivative.
#' @param kernel_type String, the kernel.
#' @param regularize Boolean.
#' @return bandwidth sequence
#' @keywords internal
bw_rot = function(y_data, x_data, y_grid, x, p, q, mu, nu, kernel_type, regularize){
  sd_y = stats::sd(y_data)
  sd_x = apply(x_data, 2, stats::sd)
  mx = apply(x_data, 2, mean)
  my = mean(y_data)
  d = ncol(x_data)
  n = length(y_data)
  ng = length(y_grid)
  data = cbind(y_data, x_data)
  if (d==1){
    # bias estimate, no rate added, DGP constant
    bx = 1.06*n^(-1/5)*sd_x
    bias_dgp = matrix(NA, ncol=3, nrow=ng)
    lower_x = min(x_data)-x
    upper_x = max(x_data)-x
    for (j in 1:ng) {
      lower_y = min(y_data) - y_grid[j]
      upper_y = max(y_data) - y_grid[j]
      # using equation from the matrix cookbook
      bias_dgp[j, 1] = normal_dgps(y_grid[j], mu, my, sd_y) * normal_dgps(x, 2, mx, sd_x)
      bias_dgp[j, 2] = normal_dgps(y_grid[j], p+1, my, sd_y) * normal_dgps(x, 0, mx, sd_x)

      Sx = solve(S_exact(lower= lower_x, upper= upper_x,
                         eval_pt=x, p=q, kernel_type=kernel_type))
      Sy = solve(S_exact(lower= lower_y, upper=upper_y,
                         eval_pt=y_grid[j], p=p, kernel_type=kernel_type))
      cx = c_exact(lower= lower_x, upper= upper_x,eval_pt=x, m=q+1, p=q, kernel_type=kernel_type)
      cy = c_exact(lower= lower_y, upper=upper_y,eval_pt=y_grid[j], m=p+1, p=p, kernel_type=kernel_type)
      Ty = T_y_exact(lower= lower_y, upper=upper_y,eval_pt=y_grid[j], p=p)
      # TODO: substitute exact Tx function
      Tx = T_x(x_data=x_data, eval_pt=x, q=q, h = bx,
               kernel_type=kernel_type)/bx
      # Tx = T_y_exact(eval_pt=x, p=q)

      bias_dgp[j, 1] = bias_dgp[j, 1] * (t(cx)%*%Sx)[nu+1]
      bias_dgp[j, 2] = bias_dgp[j, 2] * (t(cy)%*%Sy)[mu+1]
      bias_dgp[j, 3] = (bias_dgp[j, 1] + bias_dgp[j, 2])^2
    }

    # variance estimate. See Lemma 7 in the Appendix.
    v_dgp = matrix(NA, ncol=1, nrow=ng)
    if (mu > 0){
      for (j in 1:ng) {
        Sx = solve(S_exact(lower= lower_x, upper= upper_x, eval_pt=x, p=q, kernel_type=kernel_type))
        Sy = solve(S_exact(lower= lower_y, upper=upper_y,
                         eval_pt=y_grid[j], p=p, kernel_type=kernel_type))
        Ty = T_y_exact(lower= lower_y, upper=upper_y,eval_pt=y_grid[j], p=p)
        if (mu==1){
          Tx = T_x(x_data=x_data, eval_pt=x, q=q, h = bx,
                   kernel_type=kernel_type)/bx
        }else {
          Tx = T_y_exact(lower= lower_x, upper= upper_x,eval_pt = x, p=q)
        }

        v_dgp[j, 1] = stats::dnorm(y_grid[j]) * stats::pnorm(x) *
          (Sy%*%Ty%*%Sy)[mu+1, mu+1] * (Sx%*%Tx%*%Sx)[nu+1, nu+1]
      }
    }else {
      for (j in 1:ng) {
        Sx = solve(S_exact(lower= lower_x, upper= upper_x, eval_pt=x, p=q, kernel_type=kernel_type))
        Tx = T_x(x_data=x_data, eval_pt=x, q=q, h = bx,
                 kernel_type=kernel_type)/bx
        # Tx = T_y_exact(eval_pt = x, p=q)

        cdf_hat = stats::pnorm(y_grid[j]) * stats::pnorm(x)
        v_dgp [j, 1] = cdf_hat*(1-cdf_hat)
        v_dgp[j, 1] = v_dgp[j, 1] * (Sx%*%Tx%*%Sx)[nu+1, nu+1]
      }
    }

    # bandwidth
    if (mu==0){
      alpha = d + 2*min(p, q) + 2*min(mu, nu) + 2*nu + 2
    } else{
      alpha = d + 2*min(p, q) + 2*max(mu, nu) + 1
    }
    h = (abs(v_dgp/bias_dgp[, 3]))^(1/alpha)*n^(-1/alpha)

    if (regularize == TRUE){
      for (j in 1:ng){
        h[j] <- max(h[j], sort(abs(x_data-x))[min(n, 20+q+4)])
        h[j] <- max(h[j], sort(abs(y_data-y_grid[j]))[min(n, 20+p+4)])
      }
    }

  } else {

    bx = (4/(d+1))^(1/(d+4))*n^(-1/(d+4))*sd_x
    #pdf of multivariate norm
    sigma_hat = stats::cov(x_data)

    mvt_deriv = mvtnorm::dmvnorm(t(x), mean=rep(0, nrow(sigma_hat)), sigma=sigma_hat)*(solve(sigma_hat)%*%(x-mx)%*%t(x-mx)%*%solve(sigma_hat) - solve(sigma_hat))

    # theta_(3,0) expression
    theta_dd = stats::D(stats::D(expression(exp(-0.5*(y-mu)^2/sigma^2)/sqrt(2*pi*sigma^2)), "y"), "y")

    # bias estimate, no rate added, DGP constant
    bias_dgp = matrix(0, ncol=3, nrow=ng)
    for (j in 1:ng) {
      z = matrix(c(y_grid[j], x))
      mu_hat = t(stats::cor(x_data))%*%diag(d)%*%x

      # using equation from the matrix cookbook
      #theta_3,0 eval
      bias_dgp[j, 2] = eval(theta_dd, list(y=y_grid[j], mu=my, sigma=sd_y)) *
        mvtnorm::pmvnorm(upper=x, corr = stats::cor(x_data))

      Sx = solve(S_exact(eval_pt=x, p=q, kernel_type=kernel_type))
      Sy = solve(S_exact(eval_pt=y_grid[j], p=p, kernel_type=kernel_type))
      # cx = c_exact(eval_pt=x_grid[,j], m=q+1, p=q, kernel_type=kernel_type)
      cy = c_exact(eval_pt=y_grid[j], m=p+1, p=p, kernel_type=kernel_type)
      Ty = T_y_exact(eval_pt=y_grid[j], p=p)
      # need to substitute exact Tx function
      Tx = T_x(x_data=x_data, eval_pt=x, q=q, h = bx, kernel_type=kernel_type)

      mv = mvec(q+1, d)
      #theta_1,2 eval
      b1 = exp(-0.5*(y_grid[j]-my)^2/sd_y)/(sqrt(2*pi)*sd_y) * (y_grid[j]-my)/sd_y * mvt_deriv
      diag_elems = diag(b1)
      b1[lower.tri(b1, diag=TRUE)] = NA
      pxy = as.vector(b1)[!is.na(b1)]
      #NOTE: THE CODE BELOW ONLY WORKS FOR q=1 CURRENTLY
      for(i in 1:(length(mv)-d)){
        cx = c_x(x_data=x_data, eval_pt = x, m=mv[[i]], q=q, h=bx, kernel_type = kernel_type)
        bias_dgp[j, 1] = bias_dgp[j, 1] + pxy[i]*(t(cx)%*%Sx)[nu+1]
      }
      for (i in 1:d){
        cx = c_x(x_data=x_data, eval_pt = x, m=mv[[length(mv)-d+i]], q=q, h=bx, kernel_type = kernel_type)
        bias_dgp[j, 1] = bias_dgp[j, 1] + diag_elems[i]*(t(cx)%*%Sx)[nu+1]
      }
#      bias_dgp[j, 1] = bias_dgp[j, 1] * (t(cx)%*%Sx)[nu+1]
      bias_dgp[j, 2] = bias_dgp[j, 2] * (t(cy)%*%Sy)[mu+1]
      bias_dgp[j, 3] = (bias_dgp[j, 1] + bias_dgp[j, 2])^2
    }

    # variance estimate. See Lemma 7 in the Appendix.
    v_dgp = matrix(0, ncol=1, nrow=ng)
    if (mu > 0){
      for (j in 1:ng) {
        Sx = solve(S_exact(eval_pt=x, p=q, kernel_type=kernel_type))
        Sy = solve(S_exact(eval_pt=y_grid[j], p=p, kernel_type=kernel_type))
        Ty = T_y_exact(eval_pt=y_grid[j], p=p)
        # need to substitute exact Tx function
        Tx = T_x(x_data=x_data, eval_pt=x, q=q, h = bx, kernel_type=kernel_type)

        v_dgp[j, 1] = mvtnorm::dmvnorm(matrix(c(y_grid[j],x), nrow=1), mean=rep(0, d+1), sigma=stats::cov(data))^2
        v_dgp [j, 1] = v_dgp[j, 1] * (Sy%*%Ty%*%Sy)[mu+1, mu+1] * (Sx%*%Tx%*%Sx)[nu+1, nu+1]
      }
    }else {
      for (j in 1:ng) {
        Sx = solve(S_exact(eval_pt=x, p=q, kernel_type=kernel_type))
        Sy = solve(S_exact(eval_pt=y_grid[j], p=p, kernel_type=kernel_type))
        Ty = T_y_exact(eval_pt=y_grid[j], p=p)
        # need to substitute exact tx function
        Tx = T_x(x_data=x_data, eval_pt=x, q=q, h = bx, kernel_type=kernel_type)

        cdf_hat =  mvtnorm::pmvnorm(c(y_grid[j],x), mean=rep(0, d+1), sigma=stats::cov(data))
        v_dgp [j , 1] = cdf_hat*(1-cdf_hat) * mvtnorm::dmvnorm(c(y_grid[j],x), mean=rep(0, d+1), sigma=stats::cov(data))
      }
      v_dgp[, 1] = v_dgp * (Sx%*%Tx%*%Sx)[nu+1, nu+1]
    }

    h = (v_dgp/bias_dgp[, 3])^(1/6)*n^(-1/6)
    h = stats::sd(y_data)*stats::sd(x_data)*h

  }
    return(h)
}

#######################################################################################
#' IROT Bandwidth selection
#'
#' Internal Function
#' @param y_data Numeric matrix/data frame, the raw data of independent.
#' @param x_data Numeric matrix/data frame, the raw data of covariates.
#' @param y_grid Numeric vector, the evaluation points.
#' @param x Numeric, specifies the evaluation point(s) in the x-direction.
#' @param p Integer, polynomial order.
#' @param q Integer, polynomial order.
#' @param mu Integer, order of derivative.
#' @param nu Integer, order of derivative.
#' @param kernel_type String, the kernel.
#' @param regularize Boolean.
#' @return bandwidth sequence
#' @keywords internal
bw_irot = function(y_data, x_data, y_grid, x, p, q, mu, nu, kernel_type, regularize){
  sd_y = stats::sd(y_data)
  sd_x = apply(x_data, 2, stats::sd)
  mx = apply(x_data, 2, mean)
  my = mean(y_data)
  d = ncol(x_data)
  n = length(y_data)
  ng = length(y_grid)
  data = cbind(y_data, x_data)
  if (d==1){
    # bias estimate, no rate added, DGP constant
    bx = 0.5
    bias_dgp = matrix(NA, ncol=3, nrow=ng)
    for (j in 1:ng) {
      lower_x = min(x_data)-x
      lower_y = min(y_data) - y_grid[j]
      upper_x = max(x_data) - x
      upper_y = max(y_data) - y_grid[j]
      # using equation from the matrix cookbook
      bias_dgp[j, 1] = normal_dgps(y_grid[j], mu, my, sd_y) * normal_dgps(x, 2, mx, sd_x)
      bias_dgp[j, 2] = normal_dgps(y_grid[j], p+1, my, sd_y) * normal_dgps(x, 0, mx, sd_x)

      Sx = solve(S_exact(lower= lower_x, upper= upper_x,
                         eval_pt=x, p=q, kernel_type=kernel_type))
      Sy = solve(S_exact(lower= lower_y, upper=upper_y,
                         eval_pt=y_grid[j], p=p, kernel_type=kernel_type))
      cx = c_exact(lower= lower_x, upper= upper_x,eval_pt=x, m=q+1, p=q, kernel_type=kernel_type)
      cy = c_exact(lower= lower_y, upper=upper_y,eval_pt=y_grid[j], m=p+1, p=p, kernel_type=kernel_type)
      Ty = T_y_exact(lower= lower_y, upper=upper_y,eval_pt=y_grid[j], p=p)
      # TODO: substitute exact Tx function
      Tx = T_x(x_data=x_data, eval_pt=x, q=q, h = bx,
               kernel_type=kernel_type)/bx
      # Tx = T_y_exact(eval_pt=x, p=q)

      bias_dgp[j, 1] = bias_dgp[j, 1] * (t(cx)%*%Sx)[nu+1]
      bias_dgp[j, 2] = bias_dgp[j, 2] * (t(cy)%*%Sy)[mu+1]
      bias_dgp[j, 3] = (bias_dgp[j, 1] + bias_dgp[j, 2])^2
    }

    # variance estimate. See Lemma 7 in the Appendix.
    v_dgp = matrix(NA, ncol=1, nrow=ng)
    if (mu > 0){
      for (j in 1:ng) {
        Sx = solve(S_exact(lower= lower_x, upper= upper_x, eval_pt=x, p=q, kernel_type=kernel_type))
        Sy = solve(S_exact(lower= lower_y, upper=upper_y,
                         eval_pt=y_grid[j], p=p, kernel_type=kernel_type))
        Ty = T_y_exact(lower= lower_y, upper=upper_y,eval_pt=y_grid[j], p=p)
        if (mu==1){
          Tx = T_x(x_data=x_data, eval_pt=x, q=q, h = bx,
                   kernel_type=kernel_type)/bx
        }else {
          Tx = T_y_exact(lower= lower_x, upper= upper_x,eval_pt = x, p=q)
        }

        v_dgp[j, 1] = stats::dnorm(y_grid[j]) * stats::pnorm(x) *
          (Sy%*%Ty%*%Sy)[mu+1, mu+1] * (Sx%*%Tx%*%Sx)[nu+1, nu+1]
      }
    }else {
      for (j in 1:ng) {
        Sx = solve(S_exact(lower= lower_x, upper= upper_x, eval_pt=x, p=q, kernel_type=kernel_type))
        Tx = T_x(x_data=x_data, eval_pt=x, q=q, h = bx,
                 kernel_type=kernel_type)/bx

        cdf_hat = stats::pnorm(y_grid[j]) * stats::pnorm(x)
        v_dgp [j, 1] = cdf_hat*(1-cdf_hat)
        v_dgp[j, 1] = v_dgp[j, 1] * (Sx%*%Tx%*%Sx)[nu+1, nu+1]
      }
    }

    # bandwidth
    if (mu==0){
      alpha = d + 2*min(p, q) + 2*min(mu, nu) + 2*nu + 2
    } else{
      alpha = d + 2*min(p, q) + 2*max(mu, nu) + 1
    }
    h = (abs(mean(v_dgp)/(2*mean(bias_dgp[, 3]))))^(1/alpha)*n^(-1/alpha)
    h = sd_y*sd_x*h

    if (regularize == TRUE){
      h <- max(h, sort(abs(x_data-x))[min(n, 20+q+1)])
      for(j in 1:ng){
        h <- max(h, sort(abs(y_data-y_grid[j]))[min(n, 20+p+1)])
      }
    }

  } else {
    # assuming product kernel

    bx = (4/(d+1))^(1/(d+4))*n^(-1/(d+4))*sd_x
    #pdf of multivariate norm
    sigma_hat = stats::cov(x_data)

    mvt_deriv = mvtnorm::dmvnorm(t(x), mean=rep(0, nrow(sigma_hat)), sigma=sigma_hat)*(solve(sigma_hat)%*%(x-mx)%*%t(x-mx)%*%solve(sigma_hat) - solve(sigma_hat))

    # theta_(3,0) expression
    theta_dd = stats::D(stats::D(expression(exp(-0.5*(y-mu)^2/sigma^2)/sqrt(2*pi*sigma^2)), "y"), "y")

    # bias estimate, no rate added, DGP constant
    bias_dgp = matrix(0, ncol=3, nrow=ng)
    for (j in 1:ng) {
      z = matrix(c(y_grid[j], x))
      mu_hat = t(stats::cor(x_data))%*%diag(d)%*%x

      # using equation from the matrix cookbook
      #theta_3,0 eval
      bias_dgp[j, 2] = eval(theta_dd, list(y=y_grid[j], mu=my, sigma=sd_y)) *
        mvtnorm::pmvnorm(upper=x, corr = stats::cor(x_data))

      Sx = solve(S_exact(eval_pt=x, p=q, kernel_type=kernel_type))
      Sy = solve(S_exact(eval_pt=y_grid[j], p=p, kernel_type=kernel_type))
      # cx = c_exact(eval_pt=x_grid[,j], m=q+1, p=q, kernel_type=kernel_type)
      cy = c_exact(eval_pt=y_grid[j], m=p+1, p=p, kernel_type=kernel_type)
      Ty = T_y_exact(eval_pt=y_grid[j], p=p)
      # need to substitute exact Tx function
      Tx = T_x(x_data=x_data, eval_pt=x, q=q, h = bx, kernel_type=kernel_type)

      mv = mvec(q+1, d)
      #theta_1,2 eval
      b1 = exp(-0.5*(y_grid[j]-my)^2/sd_y)/(sqrt(2*pi)*sd_y) * (y_grid[j]-my)/sd_y * mvt_deriv
      diag_elems = diag(b1)
      b1[lower.tri(b1, diag=TRUE)] = NA
      pxy = as.vector(b1)[!is.na(b1)]
      #NOTE: THE CODE BELOW ONLY WORKS FOR q=1 CURRENTLY
      for(i in 1:(length(mv)-d)){
        cx = c_x(x_data=x_data, eval_pt = x, m=mv[[i]], q=q, h=bx, kernel_type = kernel_type)
        bias_dgp[j, 1] = bias_dgp[j, 1] + pxy[i]*(t(cx)%*%Sx)[nu+1]
      }
      for (i in 1:d){
        cx = c_x(x_data=x_data, eval_pt = x, m=mv[[length(mv)-d+i]], q=q, h=bx, kernel_type = kernel_type)
        bias_dgp[j, 1] = bias_dgp[j, 1] + diag_elems[i]*(t(cx)%*%Sx)[nu+1]
      }
#      bias_dgp[j, 1] = bias_dgp[j, 1] * (t(cx)%*%Sx)[nu+1]
      bias_dgp[j, 2] = bias_dgp[j, 2] * (t(cy)%*%Sy)[mu+1]
      bias_dgp[j, 3] = mvtnorm::dmvnorm(t(z), mean=rep(0, nrow(stats::cov(data))), sigma=stats::cov(data))*(bias_dgp[j, 1] + bias_dgp[j, 2])^2
    }

    # variance estimate. See Lemma 7 in the Appendix.
    v_dgp = matrix(NA, ncol=1, nrow=ng)
    if (mu > 0){
      for (j in 1:ng) {
        Sx = solve(S_exact(eval_pt=x, p=q, kernel_type=kernel_type))
        Sy = solve(S_exact(eval_pt=y_grid[j], p=p, kernel_type=kernel_type))
        Ty = T_y_exact(eval_pt=y_grid[j], p=p)
        # need to substitute exact Tx function
        Tx = T_x(x_data=x_data, eval_pt=x, q=q, h = bx, kernel_type=kernel_type)

        v_dgp[j, 1] = mvtnorm::dmvnorm(matrix(c(y_grid[j],x), nrow=1), mean=rep(0, d+1), sigma=stats::cov(data))^2
        v_dgp [j, 1] = v_dgp[j, 1] * (Sy%*%Ty%*%Sy)[mu+1, mu+1] * (Sx%*%Tx%*%Sx)[nu+1, nu+1]
      }
    }else {
      for (j in 1:ng) {
        Sx = solve(S_exact(eval_pt=x, p=q, kernel_type=kernel_type))
        Sy = solve(S_exact(eval_pt=y_grid[j], p=p, kernel_type=kernel_type))
        Ty = T_y_exact(eval_pt=y_grid[j], p=p)
        # need to substitute exact tx function
        Tx = T_x(x_data=x_data, eval_pt=x, q=q, h = bx, kernel_type=kernel_type)

        cdf_hat =  mvtnorm::pmvnorm(c(y_grid[j],x), mean=rep(0, d+1), sigma=stats::cov(data))
        v_dgp [j , 1] = cdf_hat*(1-cdf_hat) * mvtnorm::dmvnorm(c(y_grid[j],x), mean=rep(0, d+1), sigma=stats::cov(data))
        v_dgp[j, 1] = v_dgp[j, 1] * (Sx%*%Tx%*%Sx)[nu+1, nu+1]
      }
    }

    h = (sum(v_dgp)/sum(bias_dgp[, 3]))^(1/6)*n^(-1/6)
    h = stats::sd(y_data)*stats::sd(x_data)*h

  }
    return(h)
}

########################################################################################
# #' MSE Bandwidth selection
# #'
# #' Internal Function
# #' @param y_data Numeric matrix/data frame, the raw data of independent.
# #' @param x_data Numeric matrix/data frame, the raw data of covariates.
# #' @param y_grid Numeric vector, the evaluation points.
# #' @param x Numeric, specifies the evaluation point(s) in the x-direction.
# #' @param p Integer, polynomial order.
# #' @param q Integer, polynomial order.
# #' @param mu Integer, order of derivative.
# #' @param nu Integer, order of derivative.
# #' @param kernel_type String, the kernel.
# #' @return bandwidth sequence
# #' @keywords internal
#bw_mse = function(y_data, x_data, y_grid, x, p, q, mu, nu, kernel_type){
  ##centering and scaling data
  #sd_y = stats::sd(y_data)
  #sd_x = apply(x_data, 2, stats::sd)
  #mx = apply(x_data, 2, mean)
  #my = mean(y_data)
  #y_data = (y_data - my)/sd_y
  #x_data = sweep(x_data, 2, mx)/sd_x
#
  #d = ncol(x_data)
  #n = length(y_data)
  #ng = length(y_grid)
#
  #bx = (4/(d+1))^(1/(d+4))*n^(-1/(d+4))*sd_x
  ## bx = 1.06*sd_x*n^(-1/5)
  #by = 1.06*sd_y*n^(-1/5)
#
  #if (d==1){
    ## bias estimate, no rate added, DGP constant
    #e_nu = basis_vec(x, q, nu)
    #e_mu = basis_vec(0, p, mu)
#
    #bias_dgp = matrix(NA, ncol=3, nrow=ng)
    #x_scaled = as.matrix((x_data-x)/bx)
    #if(check_inv(S_x(x_scaled, q, kernel_type)/(n*bx))[1]== TRUE){
      #Sx= solve(S_x(x_scaled, q, kernel_type)/(n*bx))
    #} else{
      #singular_flag = TRUE
      #Sx= matrix(0L, nrow = length(e_nu), ncol = length(e_nu))
    #}
    #cx = c_x(x_data, x, m=q+1, q, bx, kernel_type)
    #Tx = T_x(x_data=x_data, eval_pt=x, q=q, h = bx, kernel_type=kernel_type)
    #for (j in 1:ng) {
      #y_scaled = as.matrix((y_data-y_grid[j])/by)
      ## using equation from the matrix cookbook
      #bias_dgp[j, 1] = normal_dgps(y_grid[j], mu, my, sd_y) * normal_dgps(x, 2, mx, sd_x)
      #bias_dgp[j, 2] = normal_dgps(y_grid[j], p+1, my, sd_y) * normal_dgps(x, 0, mx, sd_x)
#
      #if(check_inv(solve(S_x(y_scaled, p, kernel_type=kernel_type)/(n*by)))[1]== TRUE){
        #Sy = solve(S_x(y_scaled, p, kernel_type=kernel_type)/(n*by))
      #} else{
        #singular_flag = TRUE
        #Sy= matrix(0L, nrow = length(e_mu), ncol = length(e_mu))
      #}
      #cy = c_x(y_data, y_grid[j], m=p+1, q=p, by, kernel_type=kernel_type)
      #Ty = T_y(y_scaled, y_scaled, p, kernel_type)/(choose(n, 2)*by^2)
#
      #bias_dgp[j, 1] = bias_dgp[j, 1] * (t(cx)%*%Sx)[nu+1]
      #bias_dgp[j, 2] = bias_dgp[j, 2] * (t(cy)%*%Sy)[mu+1]
      #bias_dgp[j, 3] = (bias_dgp[j, 1] + bias_dgp[j, 2])^2
    #}
#
    ## variance estimate. See Lemma 7 in the Appendix.
    #v_dgp = matrix(NA, ncol=1, nrow=ng)
    #if (mu > 0){
      #x_scaled = as.matrix((x_data-x)/bx)
      #if(check_inv(S_x(x_scaled, q, kernel_type)/(n*bx))[1]== TRUE){
        #Sx= solve(S_x(x_scaled, q, kernel_type)/(n*bx))
      #} else{
        #singular_flag = TRUE
        #Sx= matrix(0L, nrow = length(e_nu), ncol = length(e_nu))
      #}
      #Tx = T_x(x_data=x_data, eval_pt=x, q=q, h = bx, kernel_type=kernel_type)
      #for (j in 1:ng) {
        #y_scaled = as.matrix((y_data-y_grid[j])/by)
        #if(check_inv(solve(S_x(y_scaled, p, kernel_type=kernel_type)/(n*by)))[1]== TRUE){
          #Sy = solve(S_x(y_scaled, p, kernel_type=kernel_type)/(n*by))
        #} else{
          #singular_flag = TRUE
          #Sy= matrix(0L, nrow = length(e_mu), ncol = length(e_mu))
        #}
        #Ty = T_y(y_scaled, y_scaled, p, kernel_type)/(choose(n, 2)*by^2)
#
        #v_dgp[j, 1] = stats::dnorm(y_grid[j]) * stats::pnorm(x) * (Sy%*%Ty%*%Sy)[mu+1, mu+1] * (Sx%*%Tx%*%Sx)[nu+1, nu+1]
      #}
    #}else {
      #x_scaled = as.matrix((x_data-x)/bx)
      #if(check_inv(S_x(x_scaled, q, kernel_type)/(n*bx))[1]== TRUE){
        #Sx= solve(S_x(x_scaled, q, kernel_type)/(n*bx))
      #} else{
        #singular_flag = TRUE
        #Sx= matrix(0L, nrow = length(e_nu), ncol = length(e_nu))
      #}
      #Tx = T_x(x_data=x_data, eval_pt=x, q=q, h = bx, kernel_type=kernel_type)
      #for (j in 1:ng) {
        #cdf_hat = stats::pnorm(y_grid[j]) * stats::pnorm(x)
        #v_dgp [j, 1] = cdf_hat*(1-cdf_hat)
      #}
      #v_dgp[, 1] = v_dgp * (Sx%*%Tx%*%Sx)[nu+1, nu+1]
    #}
#
    ## bandwidth
    #alpha = 2 + 2*min(p, q) + 2*max(mu, nu)
    #h = (v_dgp/bias_dgp[, 3])^(1/alpha)*n^(-1/alpha)
    #h = sd_y*sd_x*h
#
  #} else {
    ##d>1
    ## bias estimate, no rate added, DGP constant
    #e_nu = basis_vec(x, q, nu)
    #e_mu = basis_vec(0, p, mu)
#
    #bias_dgp = matrix(0, ncol=3, nrow=ng)
    #for (j in 1:ng) {
      #x_scaled = sweep(x_data, 2, x)/(bx[j]^d)
      #if(check_inv(S_x(x_scaled, q, kernel_type)/(n*bx[j]))[1]== TRUE){
        #Sx= solve(S_x(x_scaled, q, kernel_type)/(n*bx[j]))
      #} else{
        #singular_flag = TRUE
        #Sx= matrix(0L, nrow = length(e_nu), ncol = length(e_nu))
      #}
      ##TODO: Fix cx estimation
      ##estimating theta hats
      ##print(normal_dgps(y_grid[j], mu, my, sd_y))
      ##print(normal_dgps(x, 2, mx, sd_x))
      #bias_dgp[j, 1] = normal_dgps(y_grid[j], mu, my, sd_y) * prod(normal_dgps(x, 2, mx, sd_x))
      #bias_dgp[j, 2] = normal_dgps(y_grid[j], p+1, my, sd_y) * prod(normal_dgps(x, 0, mx, sd_x))
      #mv = mvec(q+1, d)
      #for(i in 1:(length(mv)-d)){
        #cx = c_x(x_data=x_data, eval_pt = x, m=mv[[i]], q=q, h=bx, kernel_type = kernel_type)
        #bias_dgp[j, 1] = bias_dgp[j, 1] + (t(cx)%*%Sx)[nu+1]
      #}
      #for (i in 1:d){
        #cx = c_x(x_data=x_data, eval_pt = x, m=mv[[length(mv)-d+i]], q=q, h=bx, kernel_type = kernel_type)
        #bias_dgp[j, 1] = bias_dgp[j, 1] + (t(cx)%*%Sx)[nu+1]
      #}
      #bias_dgp[j, 1] = bias_dgp[j, 1] * (t(cx)%*%Sx)[nu+1]
#
      #y_scaled = as.matrix((y_data-y_grid[j])/by)
      #if(check_inv(solve(S_x(y_scaled, p, kernel_type=kernel_type)/(n*by)))[1]== TRUE){
        #Sy = solve(S_x(y_scaled, p, kernel_type=kernel_type)/(n*by))
      #} else{
        #singular_flag = TRUE
        #Sy= matrix(0L, nrow = length(e_mu), ncol = length(e_mu))
      #}
      #cy = c_x(y_data, y_grid[j], m=p+1, q=p, by, kernel_type=kernel_type)
      #bias_dgp[j, 2] = bias_dgp[j, 2] * (t(cy)%*%Sy)[mu+1]
#
      #bias_dgp[j, 3] = (bias_dgp[j, 1] + bias_dgp[j, 2])^2
    #}
    ##print(bias_dgp)
#
    ## variance estimate. See Lemma 7 in the Appendix.
    #v_dgp = matrix(0, ncol=1, nrow=ng)
    #if (mu > 0){
      #for (j in 1:ng) {
        #x_scaled = as.matrix((x_data-x)/bx[j])
        #if(check_inv(S_x(x_scaled, q, kernel_type)/(n*bx[j]))[1]== TRUE){
          #Sx= solve(S_x(x_scaled, q, kernel_type)/(n*bx[j]))
        #} else{
          #singular_flag = TRUE
          #Sx= matrix(0L, nrow = length(e_nu), ncol = length(e_nu))
        #}
        #Tx = T_x(x_data=x_data, eval_pt=x, q=q, h = bx[j], kernel_type=kernel_type)
        #y_scaled = as.matrix((y_data-y_grid[j])/by)
        #if(check_inv(solve(S_x(y_scaled, p, kernel_type=kernel_type)/(n*by)))[1]== TRUE){
          #Sy = solve(S_x(y_scaled, p, kernel_type=kernel_type)/(n*by))
        #} else{
          #singular_flag = TRUE
          #Sy= matrix(0L, nrow = length(e_mu), ncol = length(e_mu))
        #}
        #Ty = T_y(y_scaled, y_scaled, p, kernel_type)/(choose(n, 2)*by^2)
#
        #v_dgp[j, 1] = stats::dnorm(y_grid[j]) * stats::pnorm(x) * (Sy%*%Ty%*%Sy)[mu+1, mu+1] * (Sx%*%Tx%*%Sx)[nu+1, nu+1]
      #}
    #}else {
      #x_scaled = as.matrix((x_data-x)/bx)
      #if(check_inv(S_x(x_scaled, q, kernel_type)/(n*bx))[1]== TRUE){
        #Sx= solve(S_x(x_scaled, q, kernel_type)/(n*bx))
      #} else{
        #singular_flag = TRUE
        #Sx= matrix(0L, nrow = length(e_nu), ncol = length(e_nu))
      #}
      #Tx = T_x(x_data=x_data, eval_pt=x, q=q, h = bx, kernel_type=kernel_type)
      #for (j in 1:ng) {
        #cdf_hat = fhat(x_data=x_data, y_data=y_data, x=x, y_grid=y_grid[j], p=3, q=1,
                       #mu=0, nu=0, h=bx, kernel_type=kernel_type)
        #v_dgp [j, 1] = cdf_hat*(1-cdf_hat)
      #}
      #v_dgp[, 1] = v_dgp * (Sx%*%Tx%*%Sx)[nu+1, nu+1]
    #}
#
    #alpha = d + 2*min(p, q) + 2*max(mu, nu) + 1
    #h = (v_dgp/bias_dgp[, 3])^(1/alpha)*n^(-1/alpha)
    #h = sd_y*stats::sd(x_data)*h
  #}
  #return(h)
#}
#
########################################################################################
# #' IMSE Bandwidth selection
# #'
# #' Internal Function
# #' @param y_data Numeric matrix/data frame, the raw data of independent.
# #' @param x_data Numeric matrix/data frame, the raw data of covariates.
# #' @param y_grid Numeric vector, the evaluation points.
# #' @param x Numeric, specifies the evaluation point(s) in the x-direction.
# #' @param p Integer, polynomial order.
# #' @param q Integer, polynomial order.
# #' @param mu Integer, order of derivative.
# #' @param nu Integer, order of derivative.
# #' @param kernel_type String, the kernel.
# #' @return bandwidth sequence
# #' @keywords internal
#bw_imse = function(y_data, x_data, y_grid, x, p, q, mu, nu, kernel_type){
  ##centering and scaling data
  #sd_y = stats::sd(y_data)
  #sd_x = apply(x_data, 2, stats::sd)
  #mx = apply(x_data, 2, mean)
  #my = mean(y_data)
  #d = ncol(x_data)
  #n = length(y_data)
  #ng = length(y_grid)
#
  #bx = (4/(d+1))^(1/(d+4))*n^(-1/(d+4))*sd_x
  ## bx = 1.06*sd_x*n^(-1/5)
  #by = 1.06*sd_y*n^(-1/5)
#
  #if (d==1){
    ## bias estimate, no rate added, DGP constant
    #e_nu = basis_vec(x, q, nu)
    #e_mu = basis_vec(0, p, mu)
#
    #bias_dgp = matrix(NA, ncol=3, nrow=ng)
    #x_scaled = as.matrix((x_data-x)/bx)
    #if(check_inv(S_x(x_scaled, q, kernel_type)/(n*bx))[1]== TRUE){
      #Sx= solve(S_x(x_scaled, q, kernel_type)/(n*bx))
    #} else{
      #singular_flag = TRUE
      #Sx= matrix(0L, nrow = length(e_nu), ncol = length(e_nu))
    #}
    #cx = c_x(x_data, x, m=q+1, q, bx, kernel_type)
    #Tx = T_x(x_data=x_data, eval_pt=x, q=q, h = bx, kernel_type=kernel_type)
    #for (j in 1:ng) {
      #y_scaled = as.matrix((y_data-y_grid[j])/by)
      ## using equation from the matrix cookbook
      #bias_dgp[j, 1] = normal_dgps(y_grid[j], mu, my, sd_y) * normal_dgps(x, 2, mx, sd_x)
      #bias_dgp[j, 2] = normal_dgps(y_grid[j], p+1, my, sd_y) * normal_dgps(x, 0, mx, sd_x)
#
      #if(check_inv(solve(S_x(y_scaled, p, kernel_type=kernel_type)/(n*by)))[1]== TRUE){
        #Sy = solve(S_x(y_scaled, p, kernel_type=kernel_type)/(n*by))
      #} else{
        #singular_flag = TRUE
        #Sy= matrix(0L, nrow = length(e_mu), ncol = length(e_mu))
      #}
      #cy = c_x(y_data, y_grid[j], m=p+1, q=p, by, kernel_type=kernel_type)
      #Ty = T_y(y_scaled, y_scaled, p, kernel_type)/(choose(n, 2)*by^2)
#
      #bias_dgp[j, 1] = bias_dgp[j, 1] * (t(cx)%*%Sx)[nu+1]
      #bias_dgp[j, 2] = bias_dgp[j, 2] * (t(cy)%*%Sy)[mu+1]
      #bias_dgp[j, 3] = (bias_dgp[j, 1] + bias_dgp[j, 2])^2
    #}
#
    ## variance estimate. See Lemma 7 in the Appendix.
    #v_dgp = matrix(NA, ncol=1, nrow=ng)
    #if (mu > 0){
      #x_scaled = as.matrix((x_data-x)/bx)
      #if(check_inv(S_x(x_scaled, q, kernel_type)/(n*bx))[1]== TRUE){
        #Sx= solve(S_x(x_scaled, q, kernel_type)/(n*bx))
      #} else{
        #singular_flag = TRUE
        #Sx= matrix(0L, nrow = length(e_nu), ncol = length(e_nu))
      #}
      #Tx = T_x(x_data=x_data, eval_pt=x, q=q, h = bx, kernel_type=kernel_type)
      #for (j in 1:ng) {
        #y_scaled = as.matrix((y_data-y_grid[j])/by)
        #if(check_inv(solve(S_x(y_scaled, p, kernel_type=kernel_type)/(n*by)))[1]== TRUE){
          #Sy = solve(S_x(y_scaled, p, kernel_type=kernel_type)/(n*by))
        #} else{
          #singular_flag = TRUE
          #Sy= matrix(0L, nrow = length(e_mu), ncol = length(e_mu))
        #}
        #Ty = T_y(y_scaled, y_scaled, p, kernel_type)/(choose(n, 2)*by^2)
#
        #v_dgp[j, 1] = stats::dnorm(y_grid[j]) * stats::pnorm(x) * (Sy%*%Ty%*%Sy)[mu+1, mu+1] * (Sx%*%Tx%*%Sx)[nu+1, nu+1]
      #}
    #}else {
      #x_scaled = as.matrix((x_data-x)/bx)
      #if(check_inv(S_x(x_scaled, q, kernel_type)/(n*bx))[1]== TRUE){
        #Sx= solve(S_x(x_scaled, q, kernel_type)/(n*bx))
      #} else{
        #singular_flag = TRUE
        #Sx= matrix(0L, nrow = length(e_nu), ncol = length(e_nu))
      #}
      #Tx = T_x(x_data=x_data, eval_pt=x, q=q, h = bx, kernel_type=kernel_type)
      #for (j in 1:ng) {
        #cdf_hat = stats::pnorm(y_grid[j]) * stats::pnorm(x)
        #v_dgp [j, 1] = cdf_hat*(1-cdf_hat)
      #}
      #v_dgp[, 1] = v_dgp * (Sx%*%Tx%*%Sx)[nu+1, nu+1]
    #}
#
    ## bandwidth
    #alpha = 2 + 2*min(p, q) + 2*max(mu, nu)
    #h = (sum(v_dgp)/sum(bias_dgp[, 3]))^(1/alpha)*n^(-1/alpha)
    #h = sd_y*sd_x*h
  #} else{
#
    #alpha = 2 + 2*min(p, q) + 2*max(mu, nu)
    #h = (sum(v_dgp)/sum(bias_dgp[, 3]))^(1/alpha)*n^(-1/alpha)
    #h = sd_y*sd_x*h
  #}
   #return(h)
 #}


#######################################################################################
#' S matrix
#' Internal Function
#'
#' Generate Matrix
#' @param eval_pt evaluation point
#' @param p Nonnegative integer, polynomial order.
#' @param kernel_type String, kernel type
#' @return a (p+1)-by-(p+1) matrix
#' @keywords internal
S_exact = function(lower=-1, upper=1, eval_pt, p, kernel_type){
  if(length(eval_pt)==1){
    # generating upper and lower limits of integration
    lower_lim = max(lower, -1)
    upper_lim = min(upper, 1)
    # looping over range of polynomials for which we need to evaluate integral
    poly_list = matrix(0L, nrow = 2*p+1, ncol = 1)
    if (lower_lim < upper_lim){
      for (i in 0:(2*p)){
        # evaluating integral
        poly_list[i+1, 1] = int_val(i, lower_lim, upper_lim, kernel_type)
      }
    }
    # initialize matrix of polynomials to be filled
    poly_mat = matrix(0L, nrow = (p+1), ncol = (p+1))
    for (j in 1:(p+1)){
      # fill matrix by selecting (p+1) consecutive elemnts of vector to fill each column
      poly_mat[, j] = poly_list[j:(j+p)]
    }
    # matrix with all the denominators (factorials)
    fact_mat = factorial(c(0:p))%*%t(factorial(c(0:p)))
    # completing the S_y matrix
    s_y = stats::dnorm(eval_pt[1])*poly_mat/fact_mat

  } else {
    # generating matrix of polynomial orders
    poly_list = c(0)
    for (i in 1:p){
      vec_len = sum(basis_vec(eval_pt, p, i))
      poly_list = append(poly_list, rep(i, vec_len), after = length(poly_list))
    }
    poly_mat = matrix(0L, nrow = length(poly_list), ncol = length(poly_list))
    poly_mat[,1] = poly_list
    poly_mat[1, ] = poly_list
    # filling matrix and scaling
    for (i in 2:length(poly_list)){
      for (j in 2:length(poly_list)){
        if (i==j){
          poly_mat[i,j] = int_val(l= poly_mat[1, j] + poly_mat[i, 1], a=-1,
                                b=1, kernel_type=kernel_type)/(factorial(poly_mat[1, j])*factorial(poly_mat[i, 1]))
        } else{
          poly_mat[i, j] =int_val(l= poly_mat[1, j] + poly_mat[i, 1], a=-1,
                                b=1, kernel_type=kernel_type)^2/(factorial(poly_mat[1, j])*factorial(poly_mat[i, 1]))
        }
      }
    }
  # matrix with all the denominators (factorials)
  #fact_mat = factorial(c(0:p))%*%t(factorial(c(0:p)))
    s_y = poly_mat
  }
  # return the S matrix
  return(s_y)
}

#' @title T matrix
#' Internal Function
#' Generate Matrix
#' @param eval_pt evaluation point
#' @param p Nonnegative integer, polynomial order.
#' @param kernel_type String, kernel type
#' @return a (p+1)-by-(p+1) matrix
#' @keywords internal
T_y_exact = function(lower=-1, upper=1, eval_pt, p, kernel_type = "uniform"){
  # ##########################################
  # TODO: need to adapt for other kernels
  # ##########################################
  # initialize empty array for filling with integrated values
  v = matrix(0L, nrow = (p+1), ncol = (p+1))
  # limits of integration
  a = max(lower, -1)
  b = min(upper, 1)
  # looping over range of polynomials for which we need to integrate
  if (a < b){
    for (i in 0:p){
      for (j in 0:p){
        # integral evaulation
        if (kernel_type == "uniform"){
          # computing all numerators and denominators
          num1 = -(b^(i+j+3)-a^(i+j+3))
          denom1 = (i+1)*(i+2)*(i+j+3)
          num2 = -(a^(i+2)*(b^(j+1)-a^(j+1)))
          denom2 = ((i+2)*(j+1))
          num3 = b^(i+1)*(b^(j+2)-a^(j+2))
          denom3 = (i+2)*(j+2)
          # computing final integral value
          v[i+1,j+1] = (num1/denom1 + num2/denom2 + num3/denom3)/4
        } else if(kernel_type == "triangular"){

        } else if(kernel_type == "epanechnikov"){

        }
      }
    }
  }
  # matrix with all the denominators (factorials)
  fact_mat = factorial(c(0:p))%*%t(factorial(c(0:p)))
  # completing the S_y matrix
  t_y =v/fact_mat
  return(t_y)
}

#' @title C matrix
#' Internal Function
#' @param eval_pt evaluation point.
#' @param m derivative order.
#' @param p Nonnegative integer, polynomial order.
#' @param kernel_type String, kernel type
#' @return a (p+1)-by-1 matrix
#' @keywords internal
c_exact = function(lower=-1, upper=-1, eval_pt, m, p, kernel_type){
  if(length(eval_pt)==1){
    # initialize empty array for filling with integrated values
    v = matrix(0L, nrow = p+1, ncol = 1)
    # setting limits of integration
    lower_lim = max(lower, -1)
    upper_lim = min(upper, 1)
    if (lower_lim < upper_lim){
      for (i in 0:p){
        # evaluted integrals
        v[i+1, 1] = int_val(i+m, lower_lim, upper_lim, kernel_type)/(factorial(m))
      }
    }
  } else {
    v = sum(basis_vec(eval_pt, p, i))
  }
return(v)
}

dmvnorm_deriv1 = function(X, mu=rep(0,ncol(X)), sigma=diag(ncol(X))) {
  fn = function(x) -1 * c((1/sqrt(det(2*pi*sigma))) * exp(-0.5*t(x-mu)%*%solve(sigma)%*%(x-mu))) * solve(sigma,(x-mu))
  out = t(apply(X,1,fn))
  return(out)
}
# Internal Function
#
# Generate Matrix
# @param eval_pt evaluation point
# @param p Nonnegative integer, polynomial order.
# @param kernel_type String, kernel type
# @return a (p+1)-by-(p+1) matrix
# @keywords internal
# S_x_exact = function(eval_pt, p, kernel_type){
#   # generating upper and lower limits of integration
#   lower_lim = -1
#   upper_lim = 1
#   # looping over range of polynomials for which we need to evaluate integral
#   poly_list = matrix(0L, nrow = 2*p+1, ncol = 1)
#   if (lower_lim < upper_lim){
#     for (i in 0:(2*p)){
#       # evaluating integral
#       poly_list[i+1, 1] = int_val(i, lower_lim, upper_lim, kernel_type)
#     }
#   }
#   # initialize matrix of polynomials to be filled
#   poly_mat = matrix(0L, nrow = (p+1), ncol = (p+1))
#   for (j in 1:(p+1)){
#     # fill matrix by selecting (p+1) consecutive elemnts of vector to fill each column
#     poly_mat[, j] = poly_list[j:(j+p)]
#   }
#   # matrix with all the denominators (factorials)
#   fact_mat = factorial(c(0:p))%*%t(factorial(c(0:p)))
#   # completing the S_y matrix
#   s_y = poly_mat/fact_mat
#   # return the S_y matrix
#   return(s_y)
# }
################################################################################
#' Internal function.
#'
#' Calculates density and higher order derivatives for Gaussian models.
#'
#' @param x Scalar, point of evaluation.
#' @param v Nonnegative integer, the derivative order (0 indicates cdf, 1 indicates pdf, etc.).
#' @param mean Scalar, the mean of the normal distribution.
#' @param sd Strictly positive scalar, the standard deviation of the normal distribution.
#'
#' @keywords internal

normal_dgps = function(x, v, mean, sd) {
  if (v == 0) {
    return(stats::pnorm(x, mean=mean, sd=sd))
  } else {
    temp = expression(exp(-(x-mean)^2/(2*sd^2))/sqrt(2*pi*sd^2))
    while(v > 1) {
      temp = stats::D(temp, "x")
      v = v - 1
    }
    return(eval(temp, list(x=x, mean=mean, sd=sd)))
  }
}
