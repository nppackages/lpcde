test_that("lpbwcde default output", {
  set.seed(42)
  n=1000
  x_data = matrix(rnorm(1*n, mean=0, sd=1), ncol=1)
  y_data = matrix(rnorm(n, mean=0, sd=1))

 #bw estimation
  model1 = lpbwcde(x_data=x_data, y_data=y_data, x=0, bw_type = "imse-rot")
  print(model1)
  summary(model1)
  coef(model1)
  expect_equal(model1$opt$ng, 19)
})
