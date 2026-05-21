# Script for generating table 2 in Cattaneo, Chandak, Jansson and Ma (2022)

#parameters
n = 2000
num_sims = 100
y_points = c(0.0, 0.8, 1.5)
x_points = c(0.0, 0.8, 1.5)

pw_testing = function(s){
  # simulating data
  x_data = matrix(rnorm(n, mean=0, sd=1))
  y_data = matrix(rnorm(n, mean=x_data, sd=1))

  true_dens = stats::dnorm(y_points, mean = x, sd=1)

  #density estimation
  est =  lpcde::lpcde(y_data=y_data, x_data=x_data, y_grid=y_points, x=x, bw=bw_star, rbc=TRUE)

  # extracting values of interest
  f_hat = est$Estimate[,3]
  f_hat_rbc = est$Estimate[,4]
  bias = (f_hat - true_dens)
  bias_rbc = (f_hat_rbc - true_dens)
  sd = est$Estimate[,5]
  sd_rbc = est$Estimate[,6]
  rmse = sqrt(sd^2 +bias^2)

  # coverage probabilities
  ce = 100*ifelse((abs(bias/sd)<=1.96), 1, 0)
  ce_rbc = 100*ifelse((abs(bias_rbc/sd_rbc)<=1.96), 1, 0)

  #average width
  aw_pw = 2*1.96*sd
  aw_pw_rbc = 2*1.96*sd_rbc

  results = matrix(c(bw_star, bias, sd, rmse, ce, ce_rbc, aw_pw, aw_pw_rbc),
                   ncol = 8)
}

pw_table = list()
bw_mat = matrix(c(0.48, 0.55, 0.78, 0.65, 0.6, 0.68, 1, 0.9, 0.9), nrow=3)
# running simulations
RNGkind(kind = "L'Ecuyer-CMRG")
set.seed(42)
for (i in 1:length(x_points)){
  x = x_points[i]
  bw_star= bw_mat[, i]
  object = parallel::mclapply(1:num_sims, pw_testing, mc.cores =8)
  avg = Reduce('+', object)/num_sims
  avg[, 2] = abs(avg[,2])
  colnames(avg) = c("BW", "bias", "sd", "rmse", "WBC CE", "RBC CE",
                    "WBC AL", "RBC AL")
  pw_table[[i]] = avg
}
library(xtable)
print(xtable(pw_table[[1]]), file="pw_table_int.tex")
print(xtable(pw_table[[2]]), file="pw_table_nb.tex")
print(xtable(pw_table[[3]]), file="pw_table_b.tex")

