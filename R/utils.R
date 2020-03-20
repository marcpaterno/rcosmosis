#' Transform a dataframe carrying grid-oriented data into a list of grid
#' vectors and a data matrix.
#'
#' @param d A data.frame with 3 columns. The first two columns must be vectors of grid
#'   coordinates, and the final column contains the value at that grid coordinate.
#' @return a list of length 3; the first 2 entries are vectors of the grid values in each
#'   direction, and the last entry is the 2 dimensional matrix of the data
#'   values.
#' @export
cols2vmat <- function(d) {
  checkmate::check_data_frame(d, ncols = 3)
  d <- dplyr::arrange_all(d)
  xvals <- unique(dplyr::pull(d, 1))
  yvals <- unique(dplyr::pull(d, 2))
  vals <- dplyr::pull(d)
  res <- list(xvals, yvals, t(matrix(vals, ncol = length(xvals), nrow = length(yvals))))
  names(res) <- names(d)
  res
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
  g <- expand.grid(u[[2]], u[[1]])
  zdata <- t(u[[3]])
  res <- tibble::tibble(x = g$Var2, y = g$Var1, z = as.numeric(zdata))
  names(res) <- names(u)
  res
}
