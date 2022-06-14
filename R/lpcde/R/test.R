library("lpcde")

set.seed(42)
n=1000
x_data = matrix(rnorm(2*n, mean=0, sd=1), ncol=2)
y_data = matrix(rnorm(n, mean=0, sd=1))
y_grid = stats::quantile(y_data, seq(from=0.1, to=0.9, by=0.1))
h = c(rep(0.2, 4), rep(0.4, 5))

# density estimation
model1 = lpcde::lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=matrix(c(0, 0), ncol=2), bw=h)
summary(model1)
