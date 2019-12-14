
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
