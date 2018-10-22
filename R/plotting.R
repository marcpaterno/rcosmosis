
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
#' @export
#'
emcee.convergence.plot <-
  function(x,
           panels = TRUE,
           numbers = FALSE,
           nmin = 1,
           nmax = max(x$n),
           walkers = 1 : max(x$walker))
  {
    nrange = nmax - nmin
    data <- filter(x
                   , walker %in% walkers
                   , nmin <= n
                   , n <= nmax)
    p <- ggplot(data, aes(x = omega_m
                          , y = SIGMA_8
                          , group = walker))
    p <- p +
      geom_point(aes(color = factor(walker)
                     , alpha = (n - nmin) / nrange)) +
      geom_path(aes(color = factor(walker)
                    , alpha = (n - nmin) / nrange)
                , show.legend = FALSE)
    if (numbers)
      p <- p + geom_text(aes(label=n), size = 3, check_overlap = TRUE)
    if (panels)
      p <- p + facet_wrap(vars(walker))
    p
  }
