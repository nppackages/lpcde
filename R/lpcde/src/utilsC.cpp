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
//' @param target Target value for sum.
//' @return List of combinations that add up to target value.
//' @keywords internal
// [[Rcpp::export]]
std::vector<std::vector<int>> print_all_sumC(int target) {
  vector<vector<int>> output;
  vector<int> result;
  print_all_sum_rec(target, 0, 1, output, result);
  return output;
}

