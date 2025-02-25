test_that("polynomial vector generator works", {
  expect_equal(mvec(2, 1), 2)
})

test_that("polynomial vector", {
  expect_equal(poly_base(2, 1), matrix(c(1, 2), nrow = 1))
})

test_that("multivariate polynomial vector", {
  expect_equal(poly_base(c(1, 2), 1), c(1, 1, 2))
})

test_that("basis vector", {
  expect_equal(basis_vec(5, 2, 1), matrix(c(0, 1, 0)))
})

test_that("multivariate basis vector", {
  expect_equal(basis_vec(c(1, 2), 1, 1), matrix(c(0, 1, 1)))
})

test_that("uniform kernel", {
  expect_equal(kernel_eval(0.5, "uniform"), 0.5)
})

test_that("triangular kernel", {
  expect_equal(kernel_eval(0.5, "triangular"), 0.5)
})

test_that("epanechnikov kernel", {
  expect_equal(kernel_eval(0.5, "epanechnikov"), 0.5625)
})

test_that("matrix inv fn", {
  expect_equal(check_inv(diag(3))[1], TRUE)
})
