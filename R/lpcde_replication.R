#Replication file for inline code in Cattaneo, Chandak, Jansson and Ma (2024)

#load package
library("lpcde")

#Setting up simulation
set.seed(42)
n=1000
x_data = matrix(rnorm(n, mean=0, sd=1))
y_data = matrix(rnorm(n, mean=x_data, sd=1))
y_grid = seq(from=-2, to=2, length.out=10)

#Standard density estimation
model1 = lpcde::lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, bw=1, rbc = TRUE)
summary(model1)
#Bandwidth selection
model2 = lpcde::lpbwcde(y_data=y_data, x_data=x_data, x=0, y_grid = y_grid)
summary(model2)

#Comparing standard and regularized estimate with true density (see Fig 1)
y_grid = seq(from=-2, to=2, length.out=20)
model1 = lpcde::lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0,
                      bw=1, rbc = TRUE)
#Regularized density estimation
model_reg = lpcde::lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, bw=1, nonneg=TRUE, normalize=TRUE)

plot(1, ylim=c(0, 0.48), ylab="density",xlab="y", xlim=c(min(y_grid), max(y_grid)) )
lines(model1$Estimate[,1], model1$Estimate[,3])
lines(model1$Estimate[,1], model_reg$Estimate[,3], col=2)
lines(model1$Estimate[,1], dnorm(y_grid), col=3)
legend('topleft',lwd=1, legend=c('standard estimate', 'regularized estimate', 'true density'), col=c(1,2,3))

#A simple plot (see Fig 2)
plot(model1, CIuniform = TRUE, rbc=TRUE, xlabel="y") + ggplot2::theme(legend.position = "none")

