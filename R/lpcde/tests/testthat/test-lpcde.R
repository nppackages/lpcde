test_that("error checking", {

})

test_that("lpcde default output", {
  set.seed(42)
  n=100
  x_data = matrix(rnorm(1*n, mean=0, sd=1), ncol=1)
  y_data = matrix(rnorm(n, mean=0, sd=1))
  y_grid = stats::quantile(y_data, seq(from=0.1, to=0.9, by=0.1))

# density estimation
  model2 = lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, bw=1.8)
  print(model2)
  confint(model2)
  coef(model2)
  vcov(model2)
  plot(model2)
  summary(model2, CIuniform=TRUE)
  confint(model2, CIuniform = TRUE)
  expect_equal(model2$opt$bw_type, "user provided")
  expect_equal(model2$opt$p, 2)
  expect_equal(model2$opt$q, 1)
  expect_equal(model2$opt$mu, 1)
  expect_equal(model2$opt$nu, 0)
  expect_equal(model2$opt$kernel, "epanechnikov")
  expect_equal(model2$opt$p_RBC, 3)
  expect_equal(model2$opt$q_RBC, 2)
})

test_that("lpcde modified output", {
  set.seed(42)
  n=100
  x_data = matrix(rnorm(1*n, mean=0, sd=1), ncol=1)
  y_data = matrix(rnorm(n, mean=0, sd=1))
  y_grid = stats::quantile(y_data, seq(from=0.1, to=0.9, by=0.1))

# density estimation
  model2 = lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, mu=0, p=3, bw_type="imse-rot")
  summary(model2)
  expect_equal(model2$opt$bw_type, "imse-rot")
  expect_equal(model2$opt$p, 3)
  expect_equal(model2$opt$q, 1)
  expect_equal(model2$opt$mu, 0)
  expect_equal(model2$opt$nu, 0)
  expect_equal(model2$opt$kernel, "epanechnikov")
  expect_equal(model2$opt$p_RBC, 4)
  expect_equal(model2$opt$q_RBC, 2)
})


test_that("error checking", {
  set.seed(42)
  n=100
  x_data = matrix(rnorm(1*n, mean=0, sd=1), ncol=1)
  y_data = matrix(rnorm(n, mean=0, sd=1))
  y_grid = stats::quantile(y_data, seq(from=0.1, to=0.9, by=0.1))

# density estimation
#  model2 = lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, p=2, q=1, p_RBC=2, q_RBC=1, bw=1.8)
#  expect_equal(model2$opt$p, model2$opt$p_RBC)
#  expect_equal(model2$opt$q, model2$opt$q_RBC)

  model2 = lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, kernel_type="triangular", bw=1.8)
  expect_equal(model2$opt$kernel, "triangular")
})

test_that("lpcde multivariate output", {
  set.seed(42)
  n=100
  x_data = matrix(rnorm(2*n, mean=0, sd=1), ncol=2)
  y_data = matrix(rnorm(n, mean=0, sd=1))
  y_grid = stats::quantile(y_data, seq(from=0.1, to=0.9, by=0.1))

# density estimation
  model2 = lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=matrix(c(0, 0), ncol=2), bw_type="imse-rot")
  summary(model2, CIuniform=TRUE)
  expect_equal(model2$opt$p, 2)
  expect_equal(model2$opt$q, 1)
  expect_equal(model2$opt$mu, 1)
  expect_equal(model2$opt$nu, 0)
  expect_equal(model2$opt$kernel, "epanechnikov")
  expect_equal(model2$opt$p_RBC, 3)
  expect_equal(model2$opt$q_RBC, 2)
})
