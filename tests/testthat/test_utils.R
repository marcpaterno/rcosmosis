library(rcosmosis)
context("CosmoSIS internal utilities")

test_that("transforming data.frame to vector/matrix form works", {
  df <- data.frame(x = rep(11:14, each = 4, rep = 2),
                   y = rep(1:2, each = 2, rep = 8),
                   z = rep(21:22, rep = 8))
  df$loglike <- with(df, -exp(-x/10 - y/20 - z/100))
  df$like <- exp(df$loglike)
  norm <- sum(df$like)
  df$like <- df$like/norm
  res <- cols2vmat(df)
  expect_is(res, "list")
  expect_equal(length(res), 4)
  expect_identical(names(res), c("x", "y", "z", "like"), info = paste(names(res), collapse = " "))
  expect_equal(res$x, 11:14)
  expect_equal(res$y, 1:2)
  expect_equal(res$z, 21:22)
  expect_is(res$like, "array")
})

test_that("transforming vector/matrix to data.frame works", {
  vm <- list(a=1:4, b=1:3, c=1:5)
  num.elements <- prod(sapply(vm, FUN=length))

})
