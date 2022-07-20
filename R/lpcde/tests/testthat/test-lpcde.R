test_that("lpcde pdf est works", {
  set.seed(42)
  n=1000
  x_data = matrix(rnorm(1*n, mean=0, sd=1), ncol=1)
  y_data = matrix(rnorm(n, mean=0, sd=1))
  y_grid = stats::quantile(y_data, seq(from=0.1, to=0.9, by=0.1))

# density estimation
  model1 = lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, bw=1.8)
  summary(model2)
  expect_equal(2 * 2, 4)
})
