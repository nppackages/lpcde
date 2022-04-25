//  utilsC.cpp
//
//  Created by Rajita Chandak on 5/25/21.
//

// [[Rcpp::depends(RcppArmadillo)]]
// #include "kcirc.hpp"
#include <iostream>
#ifdef _OPENMP
  #include <omp.h>
#endif
//#include <stdio.h>
//#include <stdlib.h>
//#include <libiomp/omp.h>
//#include "Eigen/Dense"
#include <math.h>
#include <algorithm>
#include <iterator>
#include <RcppArmadillo.h>
using namespace std;


using namespace Rcpp;
//using namespace Eigen;
using namespace arma;

//sum.cpp

// [[Rcpp::export]]
double rcpp_sum(arma::colvec v){
  double sum = 0;
  for(int i=0; i<v.size(); ++i){
    sum += v[i];
  }
  return(sum);
}

//used code from https://www.educative.io/m/find-all-sum-combinations
void print_all_sum_rec(
    int target,
    int current_sum,
    int start, vector<vector<int>>& output,
    vector<int>& result) {

  if (target == current_sum) {
    output.push_back(result);
  }

  for (int i = start; i < target; ++i) {
    int temp_sum = current_sum + i;
    if (temp_sum <= target) {
      result.push_back(i);
      print_all_sum_rec(target, temp_sum, i, output, result);
      result.pop_back();

    } else {
      return;
    }
  }
}

//' @title All Sums in C++ (Internal Function)
//' @description Function that prints all combinations of natural numbers that
//' add up to target value.
//' @param target target value for sum.
//' @return list of combinations that add up to target value.
//' @keywords internal
// [[Rcpp::export]]
std::vector<std::vector<int>> print_all_sumC(int target) {
  vector<vector<int>> output;
  vector<int> result;
  print_all_sum_rec(target, 0, 1, output, result);
  return output;
}

//kernel.cpp
// kernel function adapted for cpp, multivariate product kernel

//' @title Kernel function C++ (Internal Function)
//' @description Multivariate product kernel evaluation function. Compatible with uniform,
//' triangular and epanechnikov kernels.
//' @param x evaluation point.
//' @param kernel_type kernel function type.
//' @return kernel evaluated at \code{x}.
//' @keywords internal
// [[Rcpp::export]]
double kernel_evalC(arma::vec x, String kernel_type){
  double k = 1.0;
  int n = x.size();
  //identify if there is a value greater than 1
  // if there is, the product kernel is automatically zero
  auto missing = std::find_if(x.begin(), x.end(), [](double val){return abs(val) > 1;});
  // kernel evaluation based on kernel type
  if (missing == x.end()){
    if (kernel_type == "uniform") {
      k=pow(0.5, n);

    }else if (kernel_type == "triangular"){
      for (int i=0; i<n; i++){
        k = k*(1-abs(x(i)));
      }

    }else if (kernel_type == "epanechnikov"){
      for (int i=0; i<n; i++){
        k = k*(0.75*1.0-pow(x(i),2));
      }
    }
  }else {
    k=0.0;
  }

  return (k);
}

//poly_baseC.cpp
//adapting poly_base for multivariate input in cpp

//' @title Polynomial basis expansion in C++
//' @description polynomial basis function for univariate and multivariate input.
//' @param x a number or vector.
//' @param p order of polynomial expansion.
//' @return vector of polynomial expansion of \code{x} up to order \code{p}.
//' @keywords internal
//[[Rcpp::export]]
arma::vec poly_baseC(arma::vec x, int p){
  Environment pkg = Environment::namespace_env("lpcde");
  Function f = pkg["poly_base"];
  return as<arma::vec>(f(Rcpp::Named("x", x), Rcpp::Named("p", p)));
}


//basis_vecC.cpp
//adapting basis_vec for multivariate input in cpp
//' @title unit basis vector in C++
//' @description unit vector basis function.
//' @param x a number or vector.
//' @param p a number.
//' @param mu a number.
//' @return unit vector of appropriate length with ones corresponding to entries of order \code{mu}.
//' @keywords internal
//[[Rcpp::export]]
arma::vec basis_vecC(arma::vec x, int p, int mu){
  Environment pkg = Environment::namespace_env("lpcde");
  Function f = pkg["basis_vec"];
  return as<arma::vec>(f(Rcpp::Named("x", x), Rcpp::Named("p", p), Rcpp::Named("mu", mu)));
}

// //SymatC.cpp
// //' @title Sy matrix in C++
// //' @description Sy matrix pull from R.
// //' @param y_data raw data vector.
// //' @param y evaluation point.
// //' @param p a number.
// //' @param bw banwidth.
// //' @param kernel_type type of kernel used
// //' @return Sy matrix.
// //' @keywords internal
// //[[Rcpp::export]]
// arma::mat SymatC(arma::vec y_data, double y, int p, double bw, String kernel_type){
//   Environment pkg = Environment::namespace_env("lpcde");
//   Function f = pkg["S_y"];
//   return as<arma::mat>(f(Rcpp::Named("y_data", y_data), Rcpp::Named("eval_pt", y),
//                          Rcpp::Named("p", p), Rcpp::Named("h", bw),
//                          Rcpp::Named("kernel_type", kernel_type)));
// }

//SxmatC.cpp
//' @title S matrix in C++ (Internal Function)
//' @description Sx matrix pull from R.
//' @param x_data raw data vector or matrix.
//' @param x evaluation point.
//' @param q a number.
//' @param bw banwidth.
//' @param kernel_type type of kernel used
//' @return Sx matrix.
//' @keywords internal
//[[Rcpp::export]]
arma::mat SxmatC(arma::mat x_data, int q, String kernel_type){
  Function f("S_x");
  return as<arma::mat>(f(Rcpp::Named("x_data", x_data),
                         Rcpp::Named("p", q),
                         Rcpp::Named("kernel_type", kernel_type)));
}

// rowSums.cpp
// [[Rcpp::export]]
arma::vec rowSums(const arma::mat & X){
  int nRows = X.n_rows;
  arma::vec out(nRows);
  for(int i = 0; i < nRows; i++){
    out(i) = sum(X.row(i));
  }
  return(out);
}

// intvalC.cpp
//' @title lp integral evaluation in C++ (Internal Function)
//' @description local polynomial integral evaluation calculation of elements
//' of S_y and middle integral (evaluating integral at end points)
//' @param l degree of polynomial being integrated.
//' @param a lower limit of integration.
//' @param b upper limit of integration.
//' @param kernel_type type of kernel function. Choose from "uniform", "triangular", "epanechnikov".
//' @return value of integral.
//' @keywords internal
//[[Rcpp::export]]
double int_valC(int l, double a, double b, String kernel_type){
  double v = 0.0;
  if (kernel_type == "uniform") {
    v = 0.5*(pow(b,(l+1))-pow(a,(l+1)))/(l+1);
  }else if(kernel_type == "triangular") {
    if ((a>=0) & (b>=0)){
      double num = pow(a,(l+1))*(-2 + a + (-1 + a)*l) + pow(b, (l+1))*(2 + l - b*(1+l));
      double denom = (l+1)*(l+2);
      v = num/denom;
    }else if((a<0) & (b>0)){
      double num1 = -pow(a, (l+1))*(2+a+l+a*l);
      double denom1 = 2+(3*l)+(l*l);
      double num2 = (2+l)*pow(b, (l+1)) - (l+1)*pow(b, (l+2));
      double denom2 = (l+1)*(l+2);
      v = num1/denom1 + num2/denom2;
    }else if ((a<0) & (b<0)){
      double num = -pow(a, (l + 1))*(2 + a + l + a*l) + pow(b, (l + 1))*(2 + l + b + b*l);
      double denom = (l+1)*(l+2);
      v = num/denom;
    }
  }else if(kernel_type == "epanechnikov") {
    double num = (l+1)*pow(a, (l+3))-(l+3)*pow(a,(l+1))+pow(b,(l+1))*(3+l-pow(b,2)*(l+1));
    double denom = (l+1)*(l+3);
    v = 0.75*num/denom;
  }
  return v;
}



// // Estimator function in C++
// //' @title Function estimator (in C++)
// //' @description Estimation cdf/pdf or derivatives in C++
// //' @param y_data response variable dataset, vector.
// //' @param x_data covariate dataset, vector or matrix.
// //' @param y_grid Numeric vector, specifies the grid of evaluation points along y-direction.
// //' @param x Numeric vector or matrix, specifies the grid of evaluation points along x-direction.
// //' @param p polynomial order for y.
// //' @param q polynomial order for covariates.
// //' @param h Numeric, bandwidth vector.
// //' @param mu degree of derivative with respect to y.
// //' @param nu degree of derivative with respect to x.
// //' @param kernel_type kernel function choice.
// //' @return conditional density estimate at all grid points.
// //' @keywords internal
// //[[Rcpp::export]]
// double fhatC(arma::mat x_data, arma::mat y_data, arma::vec x, arma::vec y_grid, int p, int q, int mu, int nu, arma::vec h, String kernel_type){
//   // initializing variables
//   int n = y_data.n_rows;
//   int d = x_data.n_cols;
//   int ng = y_grid.size();
//
//   // x basis vector
//   arma::vec e_nu = basis_vecC(x, q, nu);
//
//   // y basis vector
//   arma::vec e_mu (p+1);
//   e_mu[mu] = 1;
//
//   arma::vec f_hat (ng);
//   arma::vec nh_vec (ng);
//
//   arma::vec nUnique = unique(h);
//
//   if (nUnique.size()== 1){
//     double bw = h[0];
//
//     // localization of x
//     // arma::vec idx = which(abs(x_data-x)<=h);
//
//
//
//   }
//
//   return ng;
// }
