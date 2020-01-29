
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
non.sampling.columns <- function()
{
  c("like", "post", "loglike", "prior")
}

#' Find values at which given KDE takes specified values
#'
#' @export
#' @param kde a list containing `x`- and `y`- arrays, and a `z`-matrix, as from MASS::kde2d (see below)
#' @param levels the required probability content of the contours (see below)
#'
#' @return the values of `z` that define the minimum-area regions with the specified probability contents
#' Given a list `kde` containing `x` and `y` (the x- and y-coordinates of grid points) and `z`
#' (values proportional to a probability density at those grid points), we calculate the
#' minimum-area contours which contain the probability contents specified in `levels`. All
#' the values in `levels` must be between 0 and 1.
#'
find.contours <- function(kde, levels = c(0.6826895, 0.9544997, 0.9973002))
{
  norm <- sum(kde$z)
  kde$z <- kde$z/norm
  probs.sorted <- sort(kde$z, decreasing = TRUE)
  probs.cs <- cumsum(probs.sorted)
  # Get the indices of the first values greater than the given confidence levels.
  indices <- sapply(levels, function(x) which(probs.cs > x)[1])
  probs.sorted[indices]
}


#' convert output of kde2d into a dataframe
#'
#' @param u the kde2d result to be converted
#'
#' @return a dataframe with columns `x`, `y` and `z`.
#' @export
#'
vmat2df <- function(u)
{
  g = expand.grid(u$x, u$y)
  tibble::tibble( x = g$Var1, y = g$Var2, z = as.numeric(u$z))
}
