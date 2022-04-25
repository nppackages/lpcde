#######################################################################################
# This file contains code for simulating different types of datasets
#######################################################################################

#' @title Data generating function
#' @description Function that generates dataset of prescribed size and distribution.
#' Covariates are univariate.\cr
#' Types of dgp:\cr
#' X: uniform(xmin, xmax), Y: normal(x,1),\cr
#' X: uniform(xmin, xmax), Y: uniform(min(X), max(X)),\cr
#' X: normal(0,1), Y: normal(x,1),\cr
#' X: uniform(\{0.5, 1.5\}), Y: normal(0.5, 1)+normal(1.5, 1),
#' @param n number of datapoints.
#' @param xmin minimum value of covariate.
#' @param xmax maximum value of covariate.
#' @param dgp_type joint distribution from which samples are generated.
#' @return dataset of size n.
#' @keywords internal
dgp = function(n, xmin, xmax, dgp_type){
  # setting seed to generate replicable data sets.
  set.seed(42)
  if (dgp_type == "y_normal"){
    # x ~ Unif[-1,1]
    # y ~ N(x, 1)

    # simulating x and y
    x_data = as.matrix(stats::runif(n, min = xmin, max = xmax))
    y_data = as.matrix(stats::rnorm(n, mean = x_data, sd = 1))
    all_data = cbind(y_data, x_data)

  } else if (dgp_type == "y_uniform"){
    # Let x ~ Unif[-1,1]
    # Let y ~ Unif[xmin,xmax]

    x_data = as.matrix(stats::runif(n, min = xmin, max = xmax))
    y_data = as.matrix(stats::runif(n, min = x_data, max = x_data+1))
    all_data = cbind(y_data, x_data)

  } else if (dgp_type == "joint_normal"){
    # Let x ~ N(0, 1)
    # Let y ~ N(x, 1)

    # simulating x and y
    x_data = as.matrix(stats::rnorm(n, mean = 0, sd = 1))
    y_data = as.matrix(stats::rnorm(n, mean = x_data, sd = 1))
    all_data = cbind(y_data, x_data)
  } else if (dgp_type == "truncated_normal"){
    # Let x ~ Unif(0, 1)
    # Let y ~ N(x,1)/(1-N_x(-0.8))

    x_data = as.matrix(stats::runif(n, min = -1, max = 1))
    y_data = as.matrix(ifelse(x_data>-0.8, 1,0)*stats::rnorm(n, mean = x_data, sd = 1)/(1-stats::pnorm(-0.8, mean = x_data, sd = 1)))
    all_data = cbind(y_data, x_data)

  } else if (dgp_type == "mixed_normal"){
    # Let x ~ Unif{0.5, 1.5}
    # Let y ~ N(0.5, 1) + N(1.5, 1)

    pi = c(5/10, 5/10)
    mu = c(0.5, 1.5)
    sigma = c(1, 1)

    # random generation
    rmixnorm = function(n, pi, mu, sigma) {
      k = sample.int(length(pi), n, replace = TRUE, prob = pi)
      y = stats::rnorm(n, mu[k], sigma[k])
      x = mu[k]
      all_data = cbind(y,x)
      return(all_data)
    }
    all_data = rmixnorm(n, pi, mu, sigma)

  }
  # saving generated dataset to separate folder
  # filename = paste0("simul_data/", dgp_type,"_", n, "_data.csv")

  # dataset = utils::write.table(all_data, file=filename, row.names = F, col.names = F)

  return(all_data)
}

#' @title PDF evaluation
#' @description Function that evaluates true pdf at conditioned value based on joint distribution.
#' @param x conditioned covariate value.
#' @param x_star evaluation point.
#' @param dgp_type true joint distribution.
#' @return value between 0 and 1.
#' @keywords internal
dgp_eval = function(x, x_star, dgp_type){
  if (dgp_type == "y_normal"){
    # Let x ~ Unif[-1,1]
    # Let y ~ N(x, 1)
    y = stats::dnorm(x_star, mean=x, sd=1)

  } else if (dgp_type == "y_uniform"){
    #
    # y = 1/(xmax-xmin)
  } else if (dgp_type == "joint_normal"){
    # Let x ~ N(0, 1)
    # Let y ~ N(x, 1)
    y = stats::dnorm(x_star, mean=x, sd=1)
  }else if (dgp_type == "mixed_normal"){
    # Let y ~ N(0.5, 1) + N(1.5, 1)
    pi = c(5/10, 5/10)
    mu = c(0.5, 1.5)
    sigma = c(1, 1)
    # density
    dmixnorm = function(x, pi, mu, sigma) {
      k = length(pi)
      n = length(x)
      if (n==1){
        y = sum(vapply(1:k, function(i) pi[i] * stats::dnorm(x, mu[i], sigma[i]), numeric(n)))
      } else {
        y = rowSums(vapply(1:k, function(i) pi[i] * stats::dnorm(x, mu[i], sigma[i]), numeric(n)))
      }
      return(y)
    }
    y = dmixnorm(x_star, pi, mu, sigma)
  } else if (dgp_type == "truncated_normal"){
    # Let x ~ N(0, 1)
    # Let y ~
    y = (stats::dnorm(x_star, mean = x, sd = 1)) /(1-stats::dnorm(-0.8, mean = x, sd = 1))
  }
  return(y)
}
