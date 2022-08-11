test_that("lpbwcde default output", {
  set.seed(42)
  n=100
  x_data = matrix(rnorm(1*n, mean=0, sd=1), ncol=1)
  y_data = matrix(rnorm(n, mean=0, sd=1))
  y_grid = stats::quantile(y_data, seq(from=0.1, to=0.9, by=0.1))

 #bw estimation
  model1 = lpbwcde(x_data=x_data, y_data=y_data, x=0, bw_type = "imse-rot")
  print(model1)
  summary(model1)
  coef(model1)
  expect_equal(model1$opt$ng, 19)

  expect_error(vcov(model1), regexp="The vcov method does not support \"lpbwcde\" objects.")
  expect_error(confint(model1), regexp="The confint method does not support \"lpbwcde\" objects.")

  model1 = lpbwcde(x_data=x_data, y_data=y_data, x=0, bw_type = "mse-rot")
  expect_equal(model1$opt$bw_type, "mse-rot")
  model1 = lpbwcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, bw_type = "mse-dpi")
  expect_equal(model1$opt$bw_type, "mse-dpi")
  model1 = lpbwcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, bw_type = "imse-dpi")
  expect_equal(model1$opt$bw_type, "imse-dpi")
})

test_that("lpbwcde multivariate default output", {
  set.seed(42)
  n=100
  x_data = matrix(rnorm(2*n, mean=0, sd=1), ncol=2)
  y_data = matrix(rnorm(n, mean=0, sd=1))
  y_grid = stats::quantile(y_data, seq(from=0.1, to=0.9, by=0.1))
#bw estimation
  model1 = lpbwcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=matrix(c(0, 0), ncol=2), bw_type="imse-rot")
  summary(model1)
  expect_equal(model1$opt$bw_type, "imse-rot")
})

test_that("lpbwcde default output", {
  set.seed(42)
  n=100
  x_data = matrix(rnorm(1*n, mean=0, sd=1), ncol=1)
  y_data = matrix(rnorm(n, mean=0, sd=1))
  y_grid = stats::quantile(y_data, seq(from=0.1, to=0.9, by=0.1))

#bw estimation
  model1 = lpbwcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, mu=0, p=3, bw_type = "imse-rot")
  print(model1)
  summary(model1)
  coef(model1)
  expect_equal(model1$opt$bw_type, "imse-rot")
})
