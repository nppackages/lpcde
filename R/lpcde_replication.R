# Replication file for inline code in Cattaneo, Chandak, Jansson and Ma (2024)

# load package
library("lpcde")

# software article script
set.seed(42)
n=1000
x_data = matrix(rnorm(n, mean=0, sd=1))
y_data = matrix(rnorm(n, mean=x_data, sd=1))
y_grid = seq(from=-2, to=2, length.out=10)

# density estimation
model1 = lpcde::lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, bw=1)
summary(model1)

model_reg = lpcde::lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, bw=1, nonneg=TRUE, normalize=TRUE)
summary(model_reg)

plot(model1) + ggplot2::theme(legend.position="none")

model2 = lpcde::lpbwcde(y_data=y_data, x_data=x_data, x=0, y_grid = y_grid)
summary(model2)