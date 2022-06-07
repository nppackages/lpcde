# A simple replication file that illustrates how to use the lpcde package

# load package
library("lpcde")

# software article script
set.seed(42)
n=1000
x_data = matrix(rnorm(2*n, mean=0, sd=1), ncol=2)
#x_data = matrix(rnorm(n, mean=0, sd=1))
y_data = matrix(rnorm(n, mean=0, sd=1))
y_grid = stats::quantile(y_data, seq(from=0.1, to=0.9, by=0.1))

# density estimation
model1 = lpcde::lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=matrix(c(0, 0), ncol=1), bw=0.5)
model1 = lpcde::lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, bw=0.5)
summary(model1)

#customizing output
summary(model1, alpha=0.01)

# simple plot
model2 = lpcde::lpcde(x_data=x_data, y_data=y_data, x=0, bw_type = "mse-rot")
plot(model2) + ggplot2::theme(legend.position="none")

#derivative estimation
model3 = lpcde::lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, p=3, mu=2)
summary(model3)

#customize plot output
plot(model3, alpha=0.01) + ggplot2::theme(legend.position="none")

# bandwidth selection
y_grid = stats::quantile(y_data, seq(from=0.1, to=0.9, by=0.1))
model4 = lpcde::lpbwcde(y_data=y_data, x_data=x_data, x=0, y_grid = y_grid, bw_type = "mse-rot")
summary(model4)

# IMSE bandwidth selection
model5 = lpcde::lpbwcde(y_data=y_data, x_data=x_data, x=0, y_grid = y_grid, bw_type = "imse-rot")
summary(model5)
