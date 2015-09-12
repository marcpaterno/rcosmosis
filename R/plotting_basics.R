
# These are a few functions to help make generating 'pretty' log-log
# plots a bit easier.
# xyplot.

#' Convenience functions for drawing log axes with non-default tick
#' positions and labels.
#'


#' @param xbase The base to be used for the x-axis.
#' @param ybase The base to be used for the y-axis.
#' @param xeq \code{equispaced} option for x-axis (TRUE or FALSE).
#' @param yeq \code{equispaced} option for y-axis (TRUE or FALSE).
#' @return A list suitable for the \code{scales} argument of Lattice
#' plotting functions, e.g. \code{xyplot}.
scales.log.log <- function(xbase = 10, ybase = xbase, xeq = FALSE, yeq = xeq)
{
  list( x = list(log = xbase, equispaced = xeq)
      , y = list(log = ybase, equispaced = yeq)
      )
}

#' Return the a function that returns the function
#' \code{xscale.components.log10ticks}, for more convient use in
#' \code{xyplot}.
xscale.log <- function() latticeExtra::xscale.components.log10ticks

#' Return the a function that returns the function
#' \code{yscale.components.log10ticks}, for more convient use in
#' \code{xyplot}.
yscale.log <- function() latticeExtra::yscale.components.log10ticks

