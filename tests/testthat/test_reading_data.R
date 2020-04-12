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
  txt <- readLines(fname, n = 100L)
  nwalkers <- emcee.count.walkers(txt)
  expect_identical(nwalkers, 64L)
})

test_that("reading MCMC sampler output works", {
  fname <- system.file("extdata", "sampler-output-demo5.txt.xz", package = "rcosmosis")

  samples <- read.cosmosis.mcmc(fname, 5600L)
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
  expect_identical(nrow(samples), 25600L)
  expect_identical(names(samples), c("omega_m", "h0", "deltam", "alpha", "beta", "walker", "sample"))
  # Each walker should have he same range of sample values
  ranges <- samples %>% dplyr::group_by(walker) %>% dplyr::summarize(first = min(sample), last = max(sample))
  expect_equal(length(unique(ranges$first)), 1L)
  expect_equal(length(unique(ranges$last)), 1L)
  expect_equal(max(samples$sample), 400L)
  expect_equal(max(samples$walker), 64L)
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
  expect_equal(length(filenames), 32L)
  file.path(dirname, "chain_metro_20_*.txt.xz")
}

test_that("reading MH chains works", {

  samples <- read.metropolis.hastings(get_mh_fileglob())
  expect_identical(names(samples), c("omega_m", "sigma8_input", "concentration", "chain", "sample"))
  expect_equal(nrow(samples), 32L*1000L)
})

test_that("conversion of MH chains to mcmc.list works", {
  dirname <- system.file("extdata", "run20", package = "rcosmosis")
  fglob <- file.path(dirname, "chain_metro_20_*.txt.xz")
  samples <- read.metropolis.hastings(fglob)
  ml <- mcmc.list.from.metropolis.hastings(samples)
  expect_s3_class(ml, "mcmc.list")
  expect_equal(length(ml), max(samples$chain))
})

test_that("conversion of MH chains of unequal lengths to mcmc.list works", {
  samples <- tibble::tibble(x = runif(11),
                            sample = c(1:5, 1:6),
                            chain = c(rep(1,5), rep(2,6)))
  ml <- mcmc.list.from.metropolis.hastings(samples)
  expect_s3_class(ml, "mcmc.list")
  expect_equal(length(ml), 2L)
})

test_that("remove.burnin removes samples from each chain for MH", {
 samples <- read.metropolis.hastings(get_mh_fileglob())
 expect_equal(min(samples$sample), 1L)
 good_samples <- remove.burnin(samples, 500L)
 expect_equal(min(good_samples$sample), 501L)
})

test_that("remove.burnin removes samples from each walker for emcee", {
  fname <- system.file("extdata", "sampler-output-demo5.txt.xz", package = "rcosmosis")
  samples <- read.emcee(fname)
  expect_equal(min(samples$sample), 1)
  good_samples <- remove.burnin(samples, 123L)
  expect_equal(min(good_samples$sample), 124L)
})

test_that("reading Multinest output works", {
  #root <- rprojroot::is_r_package
  #fname <- root$find_file("inst/extdata/multinest/chain_multi_1.txt.xz")
  fname <- system.file("extdata/multinest/chain_multi_1.txt.xz", package = "rcosmosis")
  d <- read.multinest(fname)
  expect_s3_class(d, "tbl_df")
  expect_true(min(d$weight) >= .Machine$double.xmin)
  expect_equal(nrow(d), 12948L)
})

test_that("reading Multinest without trimming results works", {
  fname <- system.file("extdata/multinest", "chain_multi_1.txt.xz", package = "rcosmosis")
  expect_equal(file.exists(fname), TRUE)
  d <- read.multinest(fname, remove.small = FALSE)
  expect_s3_class(d, "tbl_df")
  expect_equal(nrow(d), 13498L)
  expect_true(min(d$weight) == 0)
})

test_that("reading a matter power directory works", {
  dirname <- system.file("extdata/demo_output_1/matter_power_nl", package = "rcosmosis")
  d <- make.matterpower.dataframe(dirname, "nl")
  expect_s3_class(d, "tbl_df")
  expect_equal(nrow(d), 80200)
})

test_that("reading matter power dataframe from missing directory returns empty dataframe", {
  d <- make.matterpower.dataframe("/dev/nonexistent nonsense")
  expect_s3_class(d, "tbl_df")
  expect_equal(nrow(d), 0)
})

test_that("cosmo.scan works", {
  dirname <- system.file("extdata/demo_output_1/matter_power_nl", package = "rcosmosis")
  x <- cosmo.scan(dirname, "k_h.txt")
  expect_vector(x, ptype = numeric(), size = 200)
})

test_that("reading theory dataframe works", {
  dirname <- system.file("extdata/demo_output_1/cmb_cl", package = "rcosmosis")
  filenames <- Sys.glob(paste(dirname, "*.txt", sep = "/"))
  nfiles <- length(filenames)
  d <- make.theory.dataframe(filenames)
  expect_s3_class(d, "data.frame")
  expect_equal(length(d), nfiles)
  expect_equal(length(d), 5)
  expect_named(d, c("bb", "ee", "ell", "te", "tt"), ignore.order = TRUE)
})

test_that("nchain on empty tibble", {
  empty <- tibble::tibble(x = c(), chain = c())
  expect_equal(nchain(empty), 0)
})

test_that("nchain works on normal tibble", {
  samples <- read.metropolis.hastings(get_mh_fileglob())
  expect_equal(nchain(samples), 32L)
})
