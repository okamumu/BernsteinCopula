
test_that("PDF returns correct value for uniform R", {
  R <- matrix(1/4, nrow = 2, ncol = 2)  # 2x2 uniform weights
  u <- 0.5
  v <- 0.5
  result <- bernstein_copula_pdf(u, v, R)
  expect_true(abs(result - 1) < 1e-6)  # should be uniform copula
})

test_that("PDF is non-negative", {
  R <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2, byrow = TRUE)
  uvals <- seq(0, 1, length.out = 5)
  vvals <- seq(0, 1, length.out = 5)
  for (u in uvals) {
    for (v in vvals) {
      expect_gte(bernstein_copula_pdf(u, v, R), 0)
    }
  }
})
