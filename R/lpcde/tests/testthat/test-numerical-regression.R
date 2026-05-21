make_numerical_baseline <- function() {
  set.seed(42)
  n <- 200
  x_data <- matrix(rnorm(n), ncol = 1)
  y_data <- matrix(rnorm(n, mean = x_data, sd = 1), ncol = 1)
  y_grid <- seq(-2, 2, by = 1)

  model1 <- lpcde::lpcde(
    x_data = x_data, y_data = y_data, y_grid = y_grid, x = 0, bw = 1
  )
  model_reg <- lpcde::lpcde(
    x_data = x_data, y_data = y_data, y_grid = y_grid, x = 0, bw = 1,
    nonneg = TRUE, normalize = TRUE
  )
  model_deriv <- lpcde::lpcde(
    x_data = x_data, y_data = y_data, y_grid = y_grid, x = 0,
    p = 3, mu = 2, bw = 1
  )
  model_tri_diag <- lpcde::lpcde(
    x_data = x_data, y_data = y_data, y_grid = y_grid, x = 0, bw = 1.1,
    kernel_type = "triangular", cov_flag = "diag"
  )
  model_unif_off <- lpcde::lpcde(
    x_data = x_data, y_data = y_data, y_grid = y_grid, x = 0, bw = 1.2,
    kernel_type = "uniform", cov_flag = "off"
  )

  yq <- stats::quantile(y_data, seq(0.2, 0.8, by = 0.2))
  bw_mse <- lpcde::lpbwcde(y_data = y_data, x_data = x_data, x = 0, y_grid = yq)
  bw_imse <- lpcde::lpbwcde(
    y_data = y_data, x_data = x_data, x = 0, y_grid = yq, bw_type = "imse-rot"
  )

  x2_data <- matrix(rnorm(2 * n), ncol = 2)
  y2_data <- matrix(rnorm(n, mean = x2_data[, 1] - 0.5 * x2_data[, 2], sd = 1), ncol = 1)
  y2_grid <- stats::quantile(y2_data, c(0.25, 0.5, 0.75))
  model_2d <- lpcde::lpcde(
    x_data = x2_data, y_data = y2_data, y_grid = y2_grid,
    x = matrix(c(0, 0), ncol = 2), bw = 1, cov_flag = "diag"
  )
  bw_2d <- lpcde::lpbwcde(
    y_data = y2_data, x_data = x2_data, x = matrix(c(0, 0), ncol = 2),
    y_grid = y2_grid, bw_type = "imse-rot"
  )

  mk_matrix <- function(case, prefix, mat) {
    col_id <- colnames(mat)
    if (is.null(col_id)) {
      col_id <- seq_len(ncol(mat))
    }
    data.frame(
      case = case,
      metric = paste0(
        prefix, "[", rep(seq_len(nrow(mat)), times = ncol(mat)), ",",
        rep(col_id, each = nrow(mat)), "]"
      ),
      value = as.vector(mat),
      stringsAsFactors = FALSE
    )
  }
  mk_value <- function(case, prefix, value) {
    if (is.null(dim(value))) {
      data.frame(
        case = case,
        metric = prefix,
        value = as.numeric(value),
        stringsAsFactors = FALSE
      )
    } else {
      mk_matrix(case, prefix, value)
    }
  }
  mk_bw <- function(case, mat) {
    data.frame(
      case = case,
      metric = paste0(
        "BW[", rep(seq_len(nrow(mat)), times = ncol(mat)), ",",
        rep(colnames(mat), each = nrow(mat)), "]"
      ),
      value = as.vector(mat),
      stringsAsFactors = FALSE
    )
  }

  rbind(
    mk_matrix("fixed_bw", "Estimate", model1$Estimate),
    mk_matrix("fixed_bw", "CovMat", model1$CovMat$CovMat),
    mk_matrix("fixed_bw", "CovMat_RBC", model1$CovMat$CovMat_RBC),
    mk_matrix("regularized_bw", "Estimate", model_reg$Estimate),
    mk_matrix("regularized_bw", "CovMat", model_reg$CovMat$CovMat),
    mk_matrix("regularized_bw", "CovMat_RBC", model_reg$CovMat$CovMat_RBC),
    mk_matrix("derivative_mu2", "Estimate", model_deriv$Estimate),
    mk_matrix("derivative_mu2", "CovMat", model_deriv$CovMat$CovMat),
    mk_matrix("derivative_mu2", "CovMat_RBC", model_deriv$CovMat$CovMat_RBC),
    mk_matrix("triangular_diag", "Estimate", model_tri_diag$Estimate),
    mk_matrix("triangular_diag", "CovMat", model_tri_diag$CovMat$CovMat),
    mk_matrix("triangular_diag", "CovMat_RBC", model_tri_diag$CovMat$CovMat_RBC),
    mk_matrix("uniform_off", "Estimate", model_unif_off$Estimate),
    mk_value("uniform_off", "CovMat", model_unif_off$CovMat$CovMat),
    mk_value("uniform_off", "CovMat_RBC", model_unif_off$CovMat$CovMat_RBC),
    mk_matrix("two_dimensional_diag", "Estimate", model_2d$Estimate),
    mk_matrix("two_dimensional_diag", "CovMat", model_2d$CovMat$CovMat),
    mk_matrix("two_dimensional_diag", "CovMat_RBC", model_2d$CovMat$CovMat_RBC),
    mk_bw("bandwidth_mse", bw_mse$BW),
    mk_bw("bandwidth_imse", bw_imse$BW),
    mk_bw("bandwidth_2d_imse", bw_2d$BW)
  )
}

test_that("fixed-seed numerical baseline is unchanged", {
  fixture_candidates <- c(
    testthat::test_path("..", "fixtures", "replication-r-smoke.csv"),
    file.path("..", "fixtures", "replication-r-smoke.csv"),
    file.path("R", "lpcde", "tests", "fixtures", "replication-r-smoke.csv")
  )
  fixture_path <- fixture_candidates[file.exists(fixture_candidates)][1]
  expect_false(is.na(fixture_path))

  fixture <- read.csv(fixture_path, stringsAsFactors = FALSE)
  fixture$value <- as.numeric(fixture$value)

  current <- make_numerical_baseline()
  fixture_key <- paste(fixture$case, fixture$metric, sep = "::")
  current_key <- paste(current$case, current$metric, sep = "::")

  expect_setequal(current_key, fixture_key)
  current <- current[match(fixture_key, current_key), ]

  expect_equal(current$value, fixture$value, tolerance = 1e-8)
})
