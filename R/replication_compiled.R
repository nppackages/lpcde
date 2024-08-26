#################################################################
## Replication script Cattaneo, Chandak, Jansson and Ma (2024) ##
#################################################################

#################################
## Replication for inline code ##
#################################
#load package
library("lpcde")

# Setting up simulation
RNGkind("L'Ecuyer-CMRG")
set.seed(30)
n=1000
x_data = matrix(rnorm(n, mean=0, sd=1))
y_data = matrix(rnorm(n, mean=x_data, sd=1))
y_grid = seq(from=-2, to=2, length.out=10)

#Standard density estimation
model1 = lpcde::lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, bw=1, rbc = TRUE)
summary(model1)

#Comparing standard and regularized estimate with true density (see Fig 1)
y_grid = seq(from=-2, to=2, length.out=20)
model1 = lpcde::lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0,
                      bw=1, rbc = TRUE)
#Regularized density estimation
model_reg = lpcde::lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, bw=1,
                        nonneg=TRUE, normalize=TRUE)

plot(1, ylim=c(0, 0.48), ylab="density",xlab="y",
     xlim=c(min(y_grid), max(y_grid)))
lines(model1$Estimate[,1], model1$Estimate[,3])
lines(model1$Estimate[,1], model_reg$Estimate[,3], col=2)
lines(model1$Estimate[,1], dnorm(y_grid), col=3)
legend('topleft',lwd=1,
       legend=c('standard estimate', 'regularized estimate', 'true density'),
       col=c(1,2,3))

#A simple plot (see Fig 2)
plot(model1, CIuniform = TRUE, rbc=TRUE, xlabel="y")

#Bandwidth selection
y_grid = seq(from=-2, to=2, length.out=10)
model2 = lpcde::lpbwcde(y_data=y_data, x_data=x_data, x=0, y_grid = y_grid)
summary(model2)

#####################################################################
## Replicating Table 2 in Cattaneo, Chandak, Jansson and Ma (2024) ##
#####################################################################

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

pw_table

#####################################################################
## Replicating Table 3 in Cattaneo, Chandak, Jansson and Ma (2024) ##
#####################################################################

#parameters
n = 2000
num_sims = 100
y = 0.0
x = 0.0
#true density function
true_dens = stats::dnorm(y, mean = x, sd=1)

hmse=0.48

#grid for bandwidth range
h_grid = c(0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5)*hmse

bw_testing = function(s){
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
RNGkind(kind = "L'Ecuyer-CMRG")
set.seed(42)
for (i in 1:length(h_grid)){
  h = h_grid[i]
  object = parallel::mclapply(1:num_sims, bw_testing, mc.cores=8)
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
rownames(table_data) = h_grid/hmse
#print table
table_data


##################################################################################
## Replication code for Section 4.2 in Cattaneo, Chandak, Jansson and Ma (2024) ##
##################################################################################

#Loading dataset
data("iris")
attach(iris)

#loading conditional density estimation packages
library("np")
library("hdrcde")
library(latex2exp)

#selecting variables
y = as.matrix(iris[, 1], ncol=1)
x = as.matrix(iris[, 3], ncol=1)
x_evals = quantile(x)

#scatter plot of data with chosen conditioning values (See Fig 3)
plot(x,y, xlab = "Petal length", ylab = "Sepal length")
abline(v=x_evals[2], col="green")
abline(v=x_evals[3], col="green")
abline(v=x_evals[4], col="green")

#estimation at first conditional value
lpcde_est_q1 = lpcde::lpcde(x_data=x, y_data=y, x=x_evals[2], bw=0.9, ng=30, nonneg = TRUE, normalize = TRUE)
hdr_est_q1 = cde(x, y, x.margin=x_evals[2], y.margin = lpcde_est_q1$Estimate[,1])
npbws_q1 = npcdensbw(x, y)
np_est_q1 = fitted(npcdens(exdat=rep(x_evals[2], 30), eydat=lpcde_est_q1$Estimate[,1], bws=npbws_q1))

#estimation at second conditional value
lpcde_est_q2 = lpcde::lpcde(x_data=x, y_data=y, x=x_evals[3], ng=30, bw=0.9, nonneg = TRUE, normalize = TRUE)
hdr_est_q2 = cde(x, y, x.margin=x_evals[3], y.margin = lpcde_est_q2$Estimate[,1])
npbws_q2 = npcdensbw(x, y)
np_est_q2 = fitted(npcdens(exdat=rep(x_evals[3], 30), eydat=lpcde_est_q2$Estimate[,1], bws=npbws_q2))

#estimation at third conditional value
lpcde_est_q3 = lpcde::lpcde(x_data=x, y_data=y, x=x_evals[4], ng=30, bw=1, nonneg = TRUE, normalize = TRUE)
hdr_est_q3 = cde(x, y, x.margin=x_evals[4], y.margin = lpcde_est_q3$Estimate[,1])
npbws_q3 = npcdensbw(x, y)
np_est_q3 = fitted(npcdens(exdat=rep(x_evals[4], 30), eydat=lpcde_est_q3$Estimate[,1], bws=npbws_q3))

#plotting all CDEs
# Figure 4(a)
plot(1, main=TeX(r'($f(y|x=1.6)$)'),  xlab="Sepal length", ylab="density", ylim=c(-0.8, 1.6), xlim=c(4.8, 6.9), cex=1.5)
lines(lpcde_est_q1$Estimate[, 1], lpcde_est_q1$Estimate[,3])
lines(hdr_est_q1$y, hdr_est_q1$z, col=2)
lines(lpcde_est_q1$Estimate[, 1], np_est_q1, col=3)
legend('topright',lwd=1, legend=c('lpcde', 'hdrcde', 'np'), col=c(1,2,3))

# Figure 5(a)
plot(1, main=TeX(r'($f(y|x=4.35)$)'),  xlab="Sepal length", ylab="density", ylim=c(-0.8, 1.6), xlim=c(4.8, 6.9))
lines(lpcde_est_q2$Estimate[, 1], lpcde_est_q2$Estimate[,3])
lines(hdr_est_q2$y, hdr_est_q2$z, col=2)
lines(lpcde_est_q2$Estimate[, 1], np_est_q2, col=3)
legend('topright',lwd=1, legend=c('lpcde', 'hdrcde', 'np'), col=c(1,2,3))

# Figure 6(a)
plot(1, main=TeX(r'($f(y|x=5.1)$)'),  xlab="Sepal length", ylab="density", ylim=c(-0.8, 1.6), xlim=c(4.8, 6.9))
lines(lpcde_est_q3$Estimate[, 1], lpcde_est_q3$Estimate[,3])
lines(hdr_est_q3$y, hdr_est_q3$z, col=2)
lines(lpcde_est_q3$Estimate[, 1], np_est_q3, col=3)
legend('topright',lwd=1, legend=c('lpcde', 'hdrcde', 'np'), col=c(1,2,3))

# Figure 4(b)
plot(lpcde_est_q1, xlabel="Sepal length", ylabel="density", title=TeX(r'($f(y|x=1.6)$ with confidence bands)'))+ ggplot2::ylim(-0.8, 1.6) + ggplot2::xlim(4.8, 6.9) + ggplot2::theme(legend.position='none', panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                                                                                                                                                                                     panel.background = ggplot2::element_blank(), text = ggplot2::element_text(size = 7))
# Figure 5(b)
plot(lpcde_est_q2,xlabel="Sepal length", ylabel="density", title=TeX(r'($f(y|x=4.35)$ with confidence bands)'))+ ggplot2::ylim(-0.8, 1.6) + ggplot2::xlim(4.8, 6.9) + ggplot2::theme(legend.position='none', panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                                                                                                                                                                                     panel.background = ggplot2::element_blank(), text = ggplot2::element_text(size = 7))
# Figure 6(b)
plot(lpcde_est_q3, xlabel="Sepal length", ylabel="density", title=TeX(r'($f(y|x=5.1)$ with confidence bands)'))+ ggplot2::ylim(-0.8, 1.6) + ggplot2::xlim(4.8, 6.9) + ggplot2::theme(legend.position='none', panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                                                                                                                                                                                     panel.background = ggplot2::element_blank(), text = ggplot2::element_text(size = 7))
