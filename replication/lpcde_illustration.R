# A simple replication file that illustrates how to use the lpcde package

# load package
library("lpcde")

# software article script
set.seed(42)
n=1000
#x_data = matrix(rnorm(2*n, mean=0, sd=1), ncol=2)
x_data = matrix(rnorm(n, mean=0, sd=1))
y_data = matrix(rnorm(n, mean=x_data, sd=1))
y_grid = seq(from=-2, to=2, by = 0.5)
# density estimation
#model1 = lpcde::lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=matrix(c(0, 0), ncol=1), bw=0.5)
model1 = lpcde::lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, bw=1)
summary(model1)
model_reg = lpcde::lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, bw=1, nonneg=TRUE, normalize=TRUE)
summary(model_reg)

set.seed(42)
n=1000
x_data = matrix(rnorm(n, mean=0, sd=1))
y_data = matrix(rnorm(n, mean=x_data, sd=1))
y_grid = seq(from=-2, to=2, by = 0.5)
# density estimation
model1 = lpcde::lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0)
model_reg = lpcde::lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, nonneg=TRUE, normalize=TRUE)
#plotting densities
#standard estimate
plot(x=y_grid, y=model1$Estimate[,3], type="l", lty=1,
     xlab="", ylab="density", ylim=c(0,0.5))
#regularized estimate
lines(x=y_grid, y=model_reg$Estimate[,3], lty=2)
#true density
lines(x=y_grid, y=dnorm(y_grid, mean = mean(y_data)), lty=1,col=2)
legend('topright',lwd=1, legend=c('f', 'normalized f', 'true f'),lty = c(1, 2, 1), col=c(1,1,2))

#customizing output
summary(model1, alpha=0.01)

# simple plot
model2 = lpcde::lpcde(x_data=x_data, y_data=y_data, x=0, y_grid=y_grid)
plot(model2) + ggplot2::theme(legend.position="none")

#derivative estimation
model3 = lpcde::lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, p=3, mu=2)
summary(model3)

#customize plot output
plot(model3, alpha=0.01) + ggplot2::theme(legend.position="none")

# bandwidth selection
print("BW selection")
y_grid = stats::quantile(y_data, seq(from=0.1, to=0.9, by=0.1))
model4 = lpcde::lpbwcde(y_data=y_data, x_data=x_data, x=0, y_grid = y_grid)
summary(model4)

# IMSE bandwidth selection
model5 = lpcde::lpbwcde(y_data=y_data, x_data=x_data, x=0, y_grid = y_grid, bw_type = "imse-rot")
summary(model5)
