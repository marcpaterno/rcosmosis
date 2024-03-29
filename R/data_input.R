#' Read data from a CosmoSIS format linearized matrix sample file, e.g. a power spectrum file.
#'
#' @param dirname The name of the directory from which we read.
#' @param filename The name of the file we will read.
#' @param quiet  Passed to \code{base::scan}.
#' @return The vector of values.
cosmo.scan <- function(dirname, filename, quiet=TRUE)
{
  scan(file.path(dirname, filename), comment.char = "#", quiet = quiet)
}

#' Create a CosmoSIS 'theory' dataframe from all the named files.
#'
#' Create a CosmoSIS 'theory' dataframe from the 'theory' files
#' generated by a run of CosmoSIS. It is assumed that each file contains
#' a single column of data, and that all the columns are of the same
#' length. This is the standard CosmoSIS output format for a variety of
#' theory calculations. Lines beginning with CosmoSIS comment character
#' "#" are ignored. The names of the columns are determined from the
#' names of the files, by dropping the file extension.
#'
#' @export
#' @param fnames A character vector containing the names of the files to be read.
#' @return A CosmoSIS 'theory' dataframe.
#
# If we change the output of CosmoSIS to make a single file with
# multiple columns, then read.table() would work directly and probably a
# bit more efficiently.
make.theory.dataframe <- function(fnames)
{
  # Scan each file, creating a vector of the right name
  columns <- lapply(fnames, function(n) scan(n, comment.char = "#", quiet = TRUE))
  result <- data.frame(columns)
  names(result) <- lapply(fnames, function(n) tools::file_path_sans_ext(basename(n)))
  result
}

#' Create a CosmoSIS matter power dataframe from all the named files.
#'
#' It is assumed the file format is the CosmoSIS format for storing power
#' spectra. Note this format is different from the format used for the storage
#' of scalar parameters.
#' @export
#' @param dirname The name of the directory containing the CosmoSIS power
#'   spectrum output.
#' @param type Any value; this value is replicated to fill the \code{type}
#'   column of the dataframe.
#' @return A CosmoSIS matter power dataframe. If the directory does not exist,
#'   the returned dataframe will be empty.
#' @examples
#' d <- make.matterpower.dataframe(here::here("inst/extdata/demo_output_1/matter_power_nl"), "nl")
#' d
make.matterpower.dataframe <- function(dirname, type)
{
  # If the directory does not exist, return an empty dataframe.
  if (!file.exists(dirname)) return(tibble::tibble())

  # Scan each file, creating a vector of the right name
  # bind the columns into a dataframe.

  z.in   <- cosmo.scan(dirname, "z.txt")
  k_h.in <- cosmo.scan(dirname, "k_h.txt")
  p_k.in <- cosmo.scan(dirname, "p_k.txt")

  # The p_k array carries the data for a matrix, with 'z' varying slowly
  # and 'k_h' varying rapidly. Thus to build the dataframe, we can rely
  # on the recycling rule to get k_h correct, but have to construct z
  # ourselves.
  dframe <- tibble::tibble(p_k = p_k.in,
                       k_h = rep(k_h.in, length(z.in)),
                       z = rep(z.in, each = length(k_h.in)))
  dframe$type = type
  dframe
}

#' Extract the parameter names from a CosmoSIS sampler output file.
#'
#' The format of the first line of CosmoSIS grid and MCMC sample output carries
#' the names of the parameters that were varied in that run of CosmoSIS. This
#' function parses that line and returns the names of the parameters.
#'
#' @param txt The text to be parsed.
#' @return A character vector containing the names of the parameters read from
#'   the file.
parse.cosmosis.parameters <- function(txt) {
  # Remove comment
  tmp <- sub("#", "", txt)
  # split on tabs
  parts <- stringr::str_split(tmp, "\t")[[1]]
  # remove leading section names
  cols <- sub("[a-zA-Z_]+--", "", parts)
  cols
  # sub("(like)|(post)", "loglike", cols, fixed = FALSE)
}

#' Create a data frame from CosmoSIS MCMC sampler output.
#'
#' Reads a file in CosmoSIS MCMC sampler output format and creates a data frame
#' from it. Each line in the file corresponds to a sample in the data frame;
#' each column in the file corresponds to a column in the data frame.
#'
#'
#' The columns in a CosmoSIS MCMC data frame are:
#' \describe{
#' \item{loglike}{log-likelihood of the sample.}
#'  \item{like}{normalized likelihood of the sample.}
#'  \item{\emph{others}}{one column per sampled variable in the MCMC
#'   output, named as in the output.}
#' }
#'
#' We expect the first line of the output to contain the names of the
#' parameters, separated by spaces, and with section names separated from
#' parameter names by a double-hyphen.
#'
#' @export
#' @param fname The name of the CosmoSIS MCMC sampler output file to be read.
#'   making the data frame
#' @param burn The length of the burn-in period; these samples are ignored in
#'   making the data frame.
#' @param drop.nonsampling (default TRUE) if TRUE, non-sampling columns are dropped
#' @return a CosmoSIS MCMC data frame
read.cosmosis.mcmc <- function(fname, burn = 0L, drop.nonsampling = TRUE)
{
  checkmate::assert_count(burn)
  checkmate::assert_file_exists(fname)
  checkmate::assert_scalar(drop.nonsampling)
  d <- utils::read.table(fname, as.is = TRUE)
  if (burn > 0L)
    d <- d[-c(1:burn), ]
  first <- readLines(fname, n = 1)
  names(d) <- parse.cosmosis.parameters(first)
  if (drop.nonsampling)
  {
    names_to_keep = setdiff(names(d), non.sampling.columns())
    d <- dplyr::select(d, dplyr::all_of(names_to_keep))
  }
  tibble::as_tibble(d)
}

#' read.emcee
#'
#' @export
#' @param fname The name of the CosmoSIS MCMC sampler output file to be read.
#'
#' @return a CosmoSIS MCMC dataframe
#'
read.emcee <- function(fname)
{
  num.walkers <- emcee.count.walkers(readLines(fname, 500L))
  x <- read.cosmosis.mcmc(fname)
  nsamples <- nrow(x)/num.walkers
  x$walker <- rep(1:num.walkers, times = nsamples)
  x$sample <- rep(1:(nsamples), each = num.walkers)
  x
}

#' Read a CosmoSIS grid sampler output file.
#'
#' Reads a file in CosmoSIS grid sampler output format and returns a list
#' describing the data.
#' @export
#' @param fname The name of the CosmoSIS MCMC sampler output file to be read.
#' @return A list of \code{n+1} components, where \code{n} is the number of
#'   coordinate axes comprising the grid.
#'
#'   \describe{ \item{\code{x}, \code{y}}{The \code{x} and \code{y} coordinates
#'   of the grid points, vectors of length \code{m} and \code{n}.}
#'   \item{\code{z}}{An \code{m} by \code{n} matrix of log-likelihoods.} }
read.cosmosis.grid <- function(fname)
{
  d <- utils::read.table(fname, as.is = TRUE)
  first <- readLines(fname, n = 1)
  names(d) <- parse.cosmosis.parameters(first)
  cols2vmat(d)
}

#' read.metropolis.hastings
#'
#' @export
#' @param fileglob The glob pattern (as used by Sys.glob) specifiying the CosmoSIS
#' Metropolis-Hastings sampler output file to be read.
#'
#' @return a CosmoSIS MCMC dataframe
#'
read.metropolis.hastings <- function(fileglob)
{
  checkmate::assert_scalar(fileglob)
  checkmate::assert_string(fileglob)
  # Determine files to be read.
  fnames = Sys.glob(fileglob)
  checkmate::assert_count(length(fnames), positive = TRUE)
  # Read all files into list of dataframes.
  tbls <- lapply(fnames, read.cosmosis.mcmc)
  # Augment each dataframe with the chain number, and the sample numbers
  chain_ids <- as.list(1:length(tbls))
  tbls_and_chain_ids <- rlist::list.zip(df = tbls, chain = chain_ids)
  tbls <- lapply(tbls_and_chain_ids,
                 function(x) {
                   ntmp <- nrow(x$df)
                   x$df |> dplyr::mutate(chain = x$chain, sample = 1:ntmp)
                   })
  # Combine into one dataframe
  dplyr::bind_rows(tbls)
}

#' read.multinest
#'
#' @param fname The name of the CosmoSIS output file to be read
#' @param remove.small If TRUE, remove very small weights
#'
#' @return a weighted CosmoSIS dataframe
#' @export
#'
read.multinest <- function(fname, remove.small = TRUE)
{
  checkmate::assert_scalar(fname)
  checkmate::assert_string(fname)
  checkmate::assert_file_exists(fname)
  checkmate::assert_scalar(remove.small)
  d <- read.cosmosis.mcmc(fname)
  d$sample <- 1:nrow(d)
  if (remove.small)
    d <- dplyr::filter(d, .data$weight > .Machine$double.xmin)
  d
}

#' emcee.count.walkers Return the number of walkers used for this EMCEE run.
#'
#' @export
#' @param txt Starting lines from the EMCEE output file
#'
#' @return The number of walkers
#'
emcee.count.walkers <- function(txt) {
  matches <- stats::na.omit(stringr::str_match(txt, "^#walkers=(\\d+)$")[,2])
  stopifnot(length(matches) == 1)
  as.integer(matches)
}

#' mcmc.list.from.emcee
#
#' @export
#' @param tbl An emcee data.frame, as created by read.emcee
#'
#' @return an mcmc.list object
#'
mcmc.list.from.emcee <- function(tbl)
{
  lst <- dplyr::group_by(tbl, .data$walker) |>
         dplyr::group_split(.keep = FALSE)
  coda::as.mcmc.list(lapply(lst, coda::as.mcmc))
}

#' mcmc.list.from.metropolis.hastings
#'
#' @param tbl A MH data.frame, as created by read.metropolis.hastings
#'
#' @return an mcmc.list object
#' @export
#' @importFrom utils head
#'
mcmc.list.from.metropolis.hastings <- function(tbl)
{
  lst <- dplyr::select(tbl, -c(sample)) |>
         dplyr::group_by(.data$chain) |>
         dplyr::group_split(.keep = FALSE)
  # Because the input dataframe might now have the same number of samples for each chain
  # (a sign that the run was terminated prematurely), truncate all to the shorted.
  min_length <- min(sapply(lst, nrow))
  lst <- lapply(lst, function(x) head(x, min_length))
  coda::as.mcmc.list(lapply(lst, coda::as.mcmc))
}

#' remove.burnin
#'
#' @export
#' @param x : a emcee or MH dataframe
#' @param n : the number of samples to remove (from each walker or chain)
#'
#' @return a copy of the input x, with n samples removed
#'
remove.burnin <- function(x, n)
{
  checkmate::assert_scalar(n)
  dplyr::filter(x, .data$sample > n)
}

nchain <- function(x)
{
  checkmate::assert_data_frame(x)
  length(unique(x$chain))
}


