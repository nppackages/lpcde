test_that("sx matrix", {
  sx = S_x(as.matrix(c(1, 2, 3)), 1, "uniform")
  sol = matrix(c(0.5, 0.5, 0.5, 0.5), nrow=2)
  expect_equal(sx, sol)
})

test_that("cx matrix", {
  sx = c_x(as.matrix(c(1, 2, 3)), 0, 2, 1, 1, "uniform")
  sol = matrix(c(0.083333333, 0.083333333), nrow=2)
  expect_equal(sx, sol)
})

