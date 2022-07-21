#######################################################################################
# This file contains code for basic utility functions used in other computations
#######################################################################################

#' @title polynomial order vector
#' @description generates list of all combinations
#' of length less than or equal to d of numbers that add up to n.
#' @param n total value of each combination
#' @param d maximum length of combinations
mvec = function(n, d){
  if (d ==1){
    mvec = n
  }else {
  pvec = print_all_sumC(n)
  for (j in 1:length(pvec)){
    if (length(pvec[[j]])<d){
      pvec[[j]] = append(pvec[[j]], rep(0, d-length(pvec[[j]])))
    }
  }
  v = c()
  for (i in 1:length(pvec)){
    v = append(v, combinat::permn(pvec[[i]]))
  }
  v = v[lengths(v)<= d]
  mvec = unique(c(v, combinat::permn(c(n, rep(0, d-1)))))
  }
  return(mvec)
}

#' @title  Polynomial basis vector expansion
#' @description Generate polynomial basis vector up to order p.
#' has multivariate functionality as described in the main paper
#' normalized by factorials in denominator.
#' NOTE: currently works only up to 4th degree polynomial expansion for multivariate \code{x}.
#' @param x a number or vector.
#' @param p a number (integer).
#' @return polynomial basis of \code{x} up to degree \code{p}.
#' @examples poly_base(x = 2, p = 5)
#' @export
poly_base = function(x, p){
  # adjust type of input
  x = as.matrix(x)
  # separate univariate and multivariate x case
  if (length(x)==1){
    # vector of ones (of length p) times x
    # change power function to go from 0 to p instead of 1 to p
    v = matrix(rep(t(x),p), ncol=p, byrow=FALSE)
    v = cbind(c(rep(1, nrow(x))), v)
    a = 0:(ncol(v)-1)

    # numerator of polynomial basis
    num = matrix(1L, nrow = length(x), ncol = p+1)
    for (i in 1:p){
      num[,i+1] = num[,i]*x
    }

    # denominator of polynomial basis
    denom = factorial(c(0:p))

    return(num/denom)
  } else {
    if (p >= 5){
      stop("Not implementable for multivariate polynomials of order greater than 4.")
    }else{
    # multivariate x case
    # generate matrix of pure exponents
    v = matrix(rep(t(x),p), ncol=p, byrow=FALSE)
    a = 0:(ncol(v)-1)
    num = matrix(1L, nrow = length(x), ncol = p+1)
    for (i in 1:p){
      num[,i+1] = num[,i]*x
    }
    num = as.matrix(num[,-1])
    # initialize polynomial vector expansion with 1 (x^0)
    polyvec = 1L
    if (p>=1){
      polyvec = c(polyvec, num[,1])
    }
    if (p>=2){
      polyvec = c(polyvec, num[,2]/matrix(factorial(2), nrow =length(num[,2]), ncol = 1 ))
      # higher order polynomial expansion terms
      for (m in 2:p){
        # get all combinations of numbers that add up to m
        factor_set = print_all_sumC(m)
        # identify number of combinations
        list_length = length(factor_set)
        # loop through combinations and evaluate for each x
        for (k in 1:list_length){
          powers = factor_set[[k]]
          freq = table(powers)
          prod_val = 0
          prod_sub_val = list()
          keepind = vector()
          for (j in 1:length(freq)){
            if (freq[[j]]>1){
              if (freq[[j]]<=length(num[,j])){
                # this does n choose k
                prod_sub_val[[j]] = utils::combn(num[,j], freq[[j]], prod)
              }else{
              }
            }else {
              # this keeps indexes for exapnsion
              keepind = c(keepind, powers[j])
            }
          }
          if (length(keepind)>0){
            keepvectors=list()
            for (l in 1:length(keepind)){
              keepvectors[[l]] = num[,keepind[l]]
            }
            keepvectors = as.data.frame(keepvectors)
            comb = purrr::cross_df(keepvectors)
            if (length(prod_sub_val)>0){
              prod_sub_val = as.data.frame(append(prod_sub_val, comb))
              prod_val = apply(purrr::cross_df(prod_sub_val), 1, prod)
            } else{
              prod_val = apply(comb, 1, prod)
            }
          } else {
            prod_val = unlist(prod_sub_val)
          }
          polyvec = c(polyvec, prod_val/matrix(factorial(m), nrow = length(prod_val), ncol = 1 ))
        }
      }
    }
    }
    return(polyvec)
  }
}

#' @title Unit basis vector
#' @description Function to generate unit basis vector according to polynomial order
#' and derivative order. This function returns unit vector that is the same size
#' as the vector returned by \code{poly_base(x, p)}.
#' @param x sample input scalar or vector.
#' @param p polynomial order.
#' @param mu derivative order.
#' @return Vector of appropriate length with ones corresponding to entries of order \code{mu}.
#' @examples basis_vec(x = 2, p = 5, mu = 1)
#' @export
basis_vec = function(x, p, mu){
  if (length(x)==1){
    e_y = matrix(0L, nrow = p+1, ncol = 1)
    e_y[mu+1] = 1
    return(e_y)
  } else {
    num = matrix(1L, nrow = length(x), ncol = p+1)
    for (i in 1:p){
      num[,i+1] = num[,i]*x
    }
    num = as.matrix(num[,-1])
    e_mu = matrix(0L, nrow = length(num[,1])+1, ncol = 1)
    if (mu == 0){
      e_mu[1]=1
    }else if (mu==1){
      e_mu = matrix(1L, nrow = length(num[,1])+1, ncol = 1)
      e_mu[1] = 0
    }
    if (p>=2){
      if (mu == 2){
        e_mu = c(e_mu,  matrix(1L, nrow = length(num[,2]), ncol = 1))
      } else {
        e_mu = c(e_mu,  matrix(0L, nrow = length(num[,2]), ncol = 1))
      }
      for (m in 2:p){
        factor_set = print_all_sumC(m)
        list_length = length(factor_set)
        for (k in 1:list_length){
          powers = factor_set[[k]]
          freq = table(powers)
          prod_val = 0
          prod_sub_val = list()
          keepind = vector()
          for (j in 1:length(freq)){
            if (freq[[j]]>1){
              if (freq[[j]]<=length(num[,j])){
                prod_sub_val[[j]] = utils::combn(num[,j], freq[[j]], prod)
              }else{
              }
            }else {
              keepind = c(keepind, powers[j])
            }
          }
          if (length(keepind)>0){
            keepvectors=list()
            for (l in 1:length(keepind)){
              keepvectors[[l]] = num[,keepind[l]]
            }
            keepvectors = as.data.frame(keepvectors)
            comb = purrr::cross_df(keepvectors)
            if (length(prod_sub_val)>0){
              prod_sub_val = as.data.frame(append(prod_sub_val, comb))
              prod_val = apply(purrr::cross_df(prod_sub_val), 1, prod)
            } else{
              prod_val = apply(comb, 1, prod)
            }
          } else {
            prod_val = unlist(prod_sub_val)
          }
          if(mu == m){
            e_mu = c(e_mu,  matrix(1L, nrow = length(prod_val), ncol = 1))
          }else{
            e_mu = c(e_mu,  matrix(0L, nrow = length(prod_val), ncol = 1))
          }
        }
      }
    }
    return(e_mu)
  }
}

#' @title lp integral (Internal Function)
#' @description local polynomial integral evaluation
#' calculation of elements of S_y and middle integral (evaluating integral at end points).
#' @param l degree of polynomial being integrated.
#' @param a lower limit of integration.
#' @param b upper limit of integration.
#' @param kernel_type type of kernel function. Choose from "uniform", "triangular", "epanechnikov".
#' @return value of integral.
#' @keywords internal
int_val = function(l, a, b, kernel_type){
  if (kernel_type == "triangular"){
    if (a>=0 && b>=0){
      num = a^(l+1)*(-2 + a + (-1 + a)*l) + b^(l+1)*(2 + l - b*(1+l))
      denom = (l+1)*(l+2)
      v = num/denom
    }else if (a<0 && b>0){
      num1 = -a^(l+1)*(2+a+l+a*l)
      denom1 = 2+3*l+l*l
      num2 = (2+l)*b^(l+1) - (l+1)*b^(l+2)
      denom2 = (l+1)*(l+2)
      v = num1/denom1 +num2/denom2
    }else if (a<0 && b<0){
      num = -a^(l + 1)*(2 + a + l + a*l) + b^(l + 1)*(2 + l + b + b*l)
      denom = (l+1)*(l+2)
      v = num/denom
    }
  }else if (kernel_type == "uniform"){
    v = 0.5*(b^(l+1)-a^(l+1))/(l+1)
  }else if (kernel_type == "epanechnikov"){
    num = (l+1)*a^(l+3)-(l+3)*a^(l+1)+b^(l+1)*(3+l-(b^2)*(l+1))
    denom = (l+1)*(l+3)
    v = 0.75*num/denom
  }
  return(v)
}

#' @title Kernel Evaluation function (Internal Function)
#' @description Function that evaluates product kernel at x based on the chosen function.
#' @param x evaluation point.
#' @param kernel_type type of kernel function. Choose from "uniform", "triangular", "epanechnikov".
#' @return kernel evaluated at \code{x}.
#' @keywords internal
kernel_eval = function(x, kernel_type){
  if (kernel_type == "uniform"){
    k = ifelse(abs(x)<=1, 0.5, 0)
  }else if (kernel_type == "triangular"){
    k = (1-abs(x))*ifelse(abs(x)<=1, 1, 0)
  }else if (kernel_type == "epanechnikov"){
    k = 0.75*(1-x^2)*ifelse(abs(x)<=1, 1, 0)
  }
  k = prod(k)
  return(k)
}

#' @title Matrix invertibility check
#' @description function to check intertibility of matrix.
#' @return TRUE if matrix is invertible.
#' @param m matrix
#' @keywords internal
check_inv  = function(m) class(try(solve(m),silent=T))=="matrix"
