
#' Transform a data.frame carrying grid-oriented data into a list of grid
#' vectors and a data matrix.
#'
#' @param d A data.frame where the first \code{n} columns are vectors of grid
#'   coordinates, and the final 2 columns contain log-likelihoods and normalized
#'   likelihoods. The grid coorindate columns should be in order from least
#'   rapidly varying to most rapidly varying.
#' @return a list of length n; the first n-1 entries are the grid values in each
#'   direction, and the last entry is the (n-1) dimensional array of the data
#'   values at each grid point.
cols2vmat <- function(d) {
  likelihoods <- d$like
  if (is.null(likelihoods)) stop(names(d))
  cols <- d[, 1:(length(d) - 2)] # Subsetting to keep all but the last two columns
  # Note that as.list does not copy its argument, it modifies it.
  cols <- as.list(cols)
  result <- lapply(cols, unique)
  ary <- array(likelihoods, dim = lapply(result, length))
  result$like <- ary
  result
}

#' Append (normalized) likelihoods to a data.frame containing log-likelihoods.
#'
#' @param d A data.frame containing a column of log-likelihoods named loglike or like.
#' @return The augmented data.frame, with a column `like` containing likelihoods
append.likelihoods <- function(d) {
  # If we have a `likes` column, test to see if it is really log(likelihood).
  # If so, rename it loglike.
  if ("like" %in% names(d))
  {
    if (any(d["like"] < 0.0)){
      # this is really log likelihood, not likelihood
      d <- dplyr::rename(d, loglike = like)
    }
  }
  likes <- exp(d$loglike)
  norm <- sum(likes)
  d$like <- likes/norm
  d
}

#' non.sampling.columns
#'
#' @return character The names of the columns in CosmoSIS output that we do *not* want to return.
#'
#' @examples
non.sampling.columns <- function()
{
  c("like", "post", "loglike", "prior")
}
