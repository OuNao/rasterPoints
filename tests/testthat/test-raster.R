test_that("data2raster executes correctly and respects output dimensions", {
  mat <- matrix(runif(1000), ncol = 2)
  usr <- c(0, 1, 0, 1)
  col_idx <- rep(1:2, 500)
  colorder <- c(1, 2)
  
  res <- data2raster(
    x = mat, col_idx = col_idx, colorder = colorder,
    usr = usr, width = 100, height = 100, cex = rep(1, 1000)
  )
  
  expect_true(is.matrix(res))
  expect_equal(dim(res), c(100, 100))
  expect_true(max(res) <= 2)
})

test_that("data2raster_density performs Min-Max normalization within bounds", {
  mat <- matrix(rnorm(2000, mean = 5), ncol = 2)
  usr <- c(range(mat[,1]), range(mat[,2]))
  
  res <- data2raster_density(
    x = mat, usr = usr, width = 100, height = 100,
    n_bins = 256, smooth = TRUE, margin_pct = 0.05
  )
  
  expect_true(is.matrix(res))
  expect_equal(dim(res), c(100, 100))
  # Values must be within 0 (empty) or 1..256 (palette range)
  expect_true(all(res >= 0 & res <= 256))
  expect_true(max(res) > 0)
})