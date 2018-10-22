library(rcosmosis)
context("Reading CosmoSIS files")

test_that("parsing parameter names from pre-fcc1161b sampler files works", {
  text <- "#cosmological_parameters--omega_m	cosmological_parameters--h0	supernova_params--deltam	supernova_params--alpha	supernova_params--beta	like"
  res <- parse.cosmosis.parameters(text)
  expect_identical(res, c("omega_m", "h0", "deltam", "alpha", "beta", "loglike"))
})

test_that("parsing parameter names from post-fcc1161b sampler files works", {
  text <- "#cosmological_parameters--omega_m	post"
  res <- parse.cosmosis.parameters(text)
  expect_identical(res, c("omega_m", "loglike"))
})

test_that("parsing number of walkers from EMCEE output works", {
  fname <- system.file("extdata", "sampler-output-demo5.txt", package = "rcosmosis")
  txt <- readLines(fname, n = 100)
  nwalkers <- emcee.count.walkers(txt)
  expect_identical(nwalkers, 64L)
})

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

test_that("reading EMCEE sampler output works", {
  fname <- system.file("extdata", "sampler-output-demo5.txt", package = "rcosmosis")

  samples <- emcee.read(fname)
  expect_is(samples, "data.frame")
  expect_identical(nrow(samples), as.integer(25600))
  expect_identical(names(samples), c("omega_m", "h0", "deltam", "alpha", "beta", "loglike", "like", "walker", "n"))
  # All the columns in samples must be numeric
  expect_equal(sum(samples$like), 1.0)
  expect_equal(max(samples$n), 400)
  expect_equal(max(samples$walker), 64)
})

test_that("reading grid sampler output works", {
  fname <- system.file("extdata", "grid-output-demo7.txt", package = "rcosmosis")

  samples <- read.cosmosis.grid(fname)
  expect_is(samples, "list")
  #expect_identical(length(samples), as.integer(3))
  #expect_identical(names(samples), c("omega_m", "sigma_8"))
})
