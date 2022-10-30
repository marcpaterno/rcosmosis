context("rcosmosis internal utilities")


test_that("cols2vmat converts 2-d grid", {
  df <- tibble::tibble(x = rep(11:14, each = 2, length.out = 8),
                       y = rep(1:2, each = 1, length.out = 8))
  df <- dplyr::mutate(df, loglike = x * y)
  res <- cols2vmat(df)
  expect_is(res, "list")
  expect_equal(length(res), 3)
  expect_named(res, c("x", "y", "loglike"))
  expect_equal(res$x, 11:14)
  expect_equal(res$y, 1:2)
  expect_is(res$loglike, "matrix")
  expect_equal(dim(res$loglike), c(4, 2))
  expect_equal(res$loglike[,1], 11:14)
  expect_equal(res$loglike[,2], 2 * (11:14))
})

test_that("vmat2df reverses cols2vmat for 2d data", {
  x <- tibble::tibble(u = rep(11:14, each = 2, length.out = 8),
                      v = rep(1:2, each = 1, length.out = 8))
  x <- dplyr::mutate(x, q = exp(-(u-12.5)**2 + (v-1.25)**2))
  vm <- cols2vmat(x)
  y <- vmat2df(vm)
  expect_is(y, "tbl_df")
  expect_equal(x, y)
})

#' Calculate fractional differences between two vectors
#'
#' @param x
#' @param y
#'
#' @return fraction difference between x and y
#'
frac_diff <- function(x, y) {
  maxes <- pmax(x, y)
  abs(x-y)/maxes
}

test_that("finding contours in unit bivarian gaussian works", {
  # Generate bivariate Gaussian samples, "x" and "y".
  set.seed(123)
  data <-
    tibble::as_tibble(mvtnorm::rmvnorm(10 * 1000,
                                       mean = c(0, 0)),
                      .name_repair = "minimal")
  names(data) <- c("x", "y")
  data <- dplyr::mutate(data, r = sqrt(x ** 2 + y ** 2))
  # Calculate and normalize the probabilities.
  kde <- MASS::kde2d(data$x, data$y, n = 51)
  kde$z <- kde$z / sum(kde$z)
  # Find the contours for 0.10, 0.30, and 0.50 percentiles.
  contours <- find.contours(kde, c(0.70, 0.30, 0.50))
  # Test that the contained fractions are correct.
  p.70 <- kde$z[kde$z >= contours[1]] |> sum
  p.30 <- kde$z[kde$z >= contours[2]] |> sum
  p.50 <- kde$z[kde$z >= contours[3]] |> sum
  expect_lte(frac_diff(p.70, 0.70), 0.01)
  expect_lte(frac_diff(p.30, 0.30), 0.01)
  expect_lte(frac_diff(p.50, 0.50), 0.01)
})

test_that("plot_density_2d can be called", {
  dirname <- system.file("extdata", "run20", package = "rcosmosis")
  fglob <- file.path(dirname, "chain_metro_20_*.txt.xz")
  samples <- read.metropolis.hastings(fglob)
  p1 <-
    dplyr::filter(samples, sample > 100) |>
    plot_density_2d(omega_m, sigma8_input)
  expect_s3_class(p1, "ggplot")
})
