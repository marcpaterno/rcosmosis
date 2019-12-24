library(rcosmosis)
context("Reading CosmoSIS files")

test_that("parsing parameter names from pre-fcc1161b sampler files works", {
  text <- "#cosmological_parameters--omega_m	cosmological_parameters--h0	supernova_params--deltam	supernova_params--alpha	supernova_params--beta	like"
  res <- parse.cosmosis.parameters(text)
  expect_identical(res, c("omega_m", "h0", "deltam", "alpha", "beta", "like"))
})

test_that("parsing parameter names from post-fcc1161b sampler files works", {
  text <- "#cosmological_parameters--omega_m	post"
  res <- parse.cosmosis.parameters(text)
  expect_identical(res, c("omega_m", "post"))
})

test_that("parsing number of walkers from EMCEE output works", {
  fname <- system.file("extdata", "sampler-output-demo5.txt.xz", package = "rcosmosis")
  txt <- readLines(fname, n = 100)
  nwalkers <- emcee.count.walkers(txt)
  expect_identical(nwalkers, 64L)
})

test_that("reading MCMC sampler output works", {
  fname <- system.file("extdata", "sampler-output-demo5.txt.xz", package = "rcosmosis")

  samples <- read.cosmosis.mcmc(fname, 5600)
  expect_s3_class(samples, "tbl_df")
  expect_identical(nrow(samples), as.integer(20*1000))
  expect_identical(names(samples), c("omega_m", "h0", "deltam", "alpha", "beta"))
  # All the columns in samples must be numeric
  expect_identical(unique(sapply(samples, class)), "numeric")
})

test_that("reading EMCEE sampler output works", {
  fname <- system.file("extdata", "sampler-output-demo5.txt.xz", package = "rcosmosis")
  samples <- read.emcee(fname)
  expect_s3_class(samples, "tbl_df")
  expect_identical(nrow(samples), as.integer(25600))
  expect_identical(names(samples), c("omega_m", "h0", "deltam", "alpha", "beta", "walker", "sample"))
  # All the columns in samples must be numeric
  expect_equal(max(samples$sample), 400)
  expect_equal(max(samples$walker), 64)
})

test_that("reading grid sampler output works", {
  fname <- system.file("extdata", "grid-output-demo7.txt", package = "rcosmosis")
  samples <- read.cosmosis.grid(fname)
  # read.cosmosis.grid returns a list, not a data.frame.
  # TODO: Reconsider this design choice.
  expect_type(samples, "list")
  expect_identical(length(samples), as.integer(3))
  expect_identical(names(samples), c("omega_m", "sigma_8", "like"))
})

test_that("conversion of EMCEE to mcmc.list works", {
  fname <- system.file("extdata", "sampler-output-demo5.txt.xz", package = "rcosmosis")
  samples <- read.emcee(fname)
  ml <- mcmc.list.from.emcee(samples)
  expect_s3_class(ml, "mcmc.list")
  expect_identical(length(ml), max(samples$walker))
})

get_mh_fileglob <- function() {
  dirname <- system.file("extdata", "run20", package = "rcosmosis")
  filenames <- dir(dirname, "chain_metro_20_.*\\.txt.xz")
  expect_equal(length(filenames), 32)
  file.path(dirname, "chain_metro_20_*.txt.xz")
}

test_that("reading MH chains works", {

  samples <- read.metropolis.hastings(get_mh_fileglob())
  expect_identical(names(samples), c("omega_m", "sigma8_input", "concentration", "chain", "sample"))
  expect_equal(nrow(samples), 32*1000)
})

test_that("conversion of MH chains to mcmc.list works", {
  dirname <- system.file("extdata", "run20", package = "rcosmosis")
  fglob <- file.path(dirname, "chain_metro_20_*.txt.xz")
  samples <- read.metropolis.hastings(fglob)
  ml <- mcmc.list.from.metropolis.hastings(samples)
  expect_s3_class(ml, "mcmc.list")
  expect_equal(length(ml), max(samples$chain))
})

test_that("remove.burnin removes samples from each chain for MH", {
 samples <- read.metropolis.hastings(get_mh_fileglob())
 expect_equal(min(samples$sample), 1)
 good_samples <- remove.burnin(samples, 500)
 expect_equal(min(good_samples$sample), 501)
})

test_that("remove.burnin removes samples from each walker for emcee", {
  fname <- system.file("extdata", "sampler-output-demo5.txt.xz", package = "rcosmosis")
  samples <- read.emcee(fname)
  expect_equal(min(samples$sample), 1)
  good_samples <- remove.burnin(samples, 123)
  expect_equal(min(good_samples$sample), 124)
})

test_that("reading Multinest output works", {
  fname <- system.file("extdata/multinest", "chain_multi_1.txt.xz", package = "rcosmosis")
  d <- read.multinest(fname)
  expect_s3_class(d, "tbl_df")
  expect_true(min(d$weight) >= .Machine$double.xmin)
  expect_equal(nrow(d), 12948)
})

test_that("reading Multinest without trimming results works", {
  fname <- system.file("extdata/multinest", "chain_multi_1.txt.xz", package = "rcosmosis")
  d <- read.multinest(fname, remove.small = FALSE)
  expect_s3_class(d, "tbl_df")
  expect_equal(nrow(d), 13498)
  expect_true(min(d$weight) == 0)
})
