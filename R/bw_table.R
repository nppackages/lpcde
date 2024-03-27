# Script for generating table 1 in Cattaneo, Chandak, Jansson and Ma (2022)

#parameters
n = 2000
num_sims = 100
y = 0.0
x = 0.0
#true density function
true_dens = stats::dnorm(y, mean = 0, sd=1)

#grid for bandwidth range
h_grid = c(0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5)*0.5

bw_testing = function(s){
  print(s)
  # simulating data
  x_data = matrix(rnorm(n, mean=0, sd=1))
  y_data = matrix(rnorm(n, mean=x_data, sd=1))

  bw_star=h
  #density estimation
  est =  lpcde::lpcde(y_data=y_data, x_data=x_data, y_grid=y, x=x, bw=bw_star, rbc=TRUE)

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

table_data = matrix(ncol=8)
# running simulations
for (i in 1:length(h_grid)){
 h = h_grid[i]
 object = parallel::mclapply(1:num_sims, bw_testing, mc.cores =8)
  int_avg = matrix(0L, nrow = num_sims, ncol = 8)
  for (l in 1:num_sims){
    int_avg[l,] = object[[l]][1,]
  }
table_data = rbind(table_data, (colSums(int_avg)/num_sims))
}
table_data = table_data[2:nrow(table_data), ]
table_data[, 2] = abs(table_data[,2])
colnames(table_data) = c("BW", "bias", "sd", "rmse", "WBC CE", "RBC CE",
                         "WBC AL", "RBC AL")
rownames(table_data) = h_grid/0.5
#print table
table_data
