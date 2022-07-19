library("lpcde")

set.seed(42)
n=1000
x_data = matrix(rnorm(1*n, mean=0, sd=1), ncol=1)
y_data = matrix(rnorm(n, mean=0, sd=1))
y_grid = stats::quantile(y_data, seq(from=0.1, to=0.9, by=0.1))

#bw estimation
model1 = lpbwcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, bw_type = "imse-rot")
summary(model1)

# density estimation
model2 = lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, bw=model1$BW[,2])
summary(model2)

#bw estimation
model1 = lpbwcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, bw_type = "mse-rot")
summary(model1)

# density estimation
model2 = lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, bw=model1$BW[,2])
summary(model2)

set.seed(42)
n=1000
x_data = matrix(rnorm(2*n, mean=0, sd=1), ncol=2)
y_data = matrix(rnorm(n, mean=0, sd=1))
y_grid = stats::quantile(y_data, seq(from=0.1, to=0.9, by=0.1))

#bw estimation
model1 = lpbwcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=matrix(c(0, 0), ncol=2), bw_type="imse-rot")
summary(model1)

# density estimation
model2 = lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=matrix(c(0, 0), ncol=2), bw=model1$BW[,2])
summary(model2)

#bw estimation
model1 = lpbwcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=matrix(c(0, 0), ncol=2), bw_type="mse-rot")
summary(model1)

# density estimation
model2 = lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=matrix(c(0, 0), ncol=2), bw=model1$BW[,2])
summary(model2)
