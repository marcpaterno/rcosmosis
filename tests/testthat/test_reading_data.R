library(rcosmosis)
context("Reading CosmoSIS files")

test_that("reading MCMC sampler output works", {
  fname <- system.file("extdata", "sampler-output-demo5.txt", package = "rcosmosis")

  samples <- read.cosmosis.mcmc(fname, 5600)
  expect_is(samples, "data.frame")
  expect_identical(nrow(samples), as.integer(20*1000))
  expect_identical(names(samples), c("omega_m", "h0", "deltam", "alpha", "beta", "loglike", "like"))
  # All the columns in samples must be numeric
  expect_identical(unique(sapply(samples, class)), "numeric")
  expect_equal(sum(samples$like), 1.0)
})

test_that("reading grid sampler otput works", {
  fname <- system.file("extdata", "grid-output-demo7.txt", package = "rcosmosis")
})
