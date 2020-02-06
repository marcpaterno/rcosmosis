#' Title emcee.convergence.plot
#'
#' @param x : a dataframe made from emcee sampler output
#' @param panels : logical, if TRUE use panels to separate walkers
#' @param numbers : logical, if TRUE print sample numbers on points
#' @param nmin : integer, first sample to plot
#' @param nmax : integer, last sample to plot
#' @param walkers : integer vector, indices of walkers to plot. If NULL, plot
#'   all walkers
#'
#' @return a ggplot object
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @export
#'
emcee.convergence.plot <-
  function(x,
           panels = TRUE,
           numbers = FALSE,
           nmin = 1,
           nmax = max(x$sample),
           walkers = 1:max(x$walker))
  {
    nrange = nmax - nmin + 1
    data <- dplyr::filter(x
                          , .data$walker %in% walkers
                          , nmin <= .data$sample
                          , .data$sample <= nmax)
    p <- ggplot(data, aes(x = .data$omega_m
                          , y = .data$SIGMA_8
                          , group = .data$walker))
    p <- p +
      geom_point(aes(
        color = factor(.data$walker)
        ,
        alpha = (.data$sample - nmin) / nrange
      )) +
      geom_path(aes(
        color = factor(.data$walker)
        ,
        alpha = (.data$sample - nmin) / nrange
      )
      , show.legend = FALSE)
    if (numbers)
      p <-
      p + geom_text(aes(label = .data$sample),
                    size = 3,
                    check_overlap = TRUE)
    if (panels)
      p <- p + facet_wrap(vars(.data$walker))
    p
  }

#' plot_density_2d
#'
#' @param data a CosmoSIS MCMC dataframe
#' @param x the variable to be plotted on the x-axis
#' @param y the variable to be plotted on the y-axis
#' @param bins the numbers of bins for both `x` and `y`
#' @param scale_x position scale function for `x`
#' @param scale_y position scale function for `y`
#' @param levels probability content to be included in each contour line
#'
#' @return a ggplot object
#'
#' @import ggplot2
#' @importFrom stats density
#' @export
#'
#' @examples
#' \dontrun{
#'   dplyr::filter(samples, sample > nburn) %>% plot_density_2d(omega_m, sigma8_input)
#' }
plot_density_2d <-
  function(data,
           x,
           y,
           bins = 50,
           scale_x = scale_x_continuous,
           scale_y = scale_y_continuous,
           levels = c(0.6826895, 0.9544997, 0.9973002))
{
  x <- enquo(x)
  y <- enquo(y)
  dxy <- dplyr::select(data, !!x, !!y)
  kde <- MASS::kde2d(dxy[[1]], dxy[[2]], n = bins)
  norm <- sum(kde$z)
  kde$z <- kde$z/norm
  breaks <- find.contours(kde, levels=levels)
  dxy <- vmat2df(kde)
  ggplot(data = data, mapping = aes(x = !!x, y = !!y)) +
    geom_hex(mapping = aes(fill = stat(log(density))), bins = bins) +
    geom_hex(mapping = aes(fill = stat(.data$count)), bins = bins) +
    scale_fill_continuous(type = "gradient") +
    geom_contour(data = dxy, mapping = aes(x, y, z=.data$z), breaks = breaks, color = "red") +
    scale_x() +
    scale_y() +
    theme_bw()
}

