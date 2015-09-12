library(DBI)
library(RSQLite)
library(lattice)
library(latticeExtra)
library(reshape2)
library(tools)

#' Create a single "distance plot". A distance plot is always a plot of
#' something as a function of reshift "z'. The plot is written to a file,
#' using the specified device. All additional arguments are passed to the
#' device function.
#' @return nothing; used for side-effect
#' @param colname As a string, the name of the column from the dataframe to be used on the y-axis.
#' @param dframe A 
make.distance.plot <- function(colname, dframe, prefix, outdir, verbose, devname, ...)
{
  device.fcn = get(devname)
  plotfile <- paste(outdir, "/", prefix, "_", colname, ".", devname, sep="")
  if (verbose) cat("Plotting", colname, "versus z into", plotfile, "\n")
  device.fcn(plotfile, ...)
  print(xyplot( as.formula(paste(colname, "~z"))
              , dframe
              , xlab="Redshift z"
              , ylab=toupper(colname)
              , grid=TRUE
              , type="l",
              , lwd = 2
              ))
  invisible(dev.off())
}

#' Create a single "CMB spectrum plot". A CMB spectrum plot is always a
#' plot of something as a function of C_ell. The plot is written to
#' a file, using the specified device. All additional arguments are
#' passed to the device function.
make.cmb.spectrum.plot <- function(colname, dframe, prefix, outdir, verbose, devname, ...)
{
  device.fcn <- get(devname)
  plotfile <- paste(outdir, "/", prefix, "_", colname, ".", devname, sep="")
  if (verbose) cat("Plotting", colname, "versus ell into", plotfile, "\n")
  device.fcn(plotfile, ...)
  print(xyplot( as.formula(paste(colname, "~ell"))
              , dframe
              , xlab = "ell"
              , ylab = paste("C_ell ", toupper(colname), "/uK^2", sep="")
              , grid=TRUE
              , type="l",
              , lwd = 2
              ))
  invisible(dev.off())
}

#' Create the CMB spectrum plot with all spectra shown together. Because
#' some of the spectra may have negative values, we take the absolute
#' value (so that our log-log plot doesn't complain).
make.cmb.grand.spectrum.plot <- function(dframe, prefix, outdir, verbose, devname, ...)
{
  device.fcn <- get(devname)
  plotfile <- paste(outdir, "/", prefix, "_grand.", devname, sep="")
  if (verbose)
    cat("Plotting all spectra versus ell into", plotfile, "\n")
  device.fcn(plotfile, ...)
  print(xyplot( abs(value)~ell
              , dframe
              , xlab="ell"
              , ylab="C_ell spectra / uK^2"
              , grid=TRUE
              , type="l",
              , lwd = 2
              , group = variable
              , auto.key = list(space="right", points=FALSE, lines=TRUE)
              , scales = scales.log.log()
              , xscale.components = xscale.log()
              , yscale.components = yscale.log()
              ))
  invisible(dev.off())
}

#' Make all the distance plots for files under the given top-level
#' directory *topdir*. Distance plots are made for data files in the
#' subdirectory "distances".
make.distance.plots <- function(topdir, verbose, prefix, outdir, devname)
{
  # data files are in "distances".
  datadir <- file.path(topdir, "distances")
  if (verbose) cat("Making distance plots from data in", datadir, "\n")

  # We don't want to read the "values.txt" file.
  excluded <- file.path(datadir, "values.txt")
  files <- Filter(function(x) {x!=excluded}, dir(datadir, full.names=TRUE))

  # Create a dataframe from all the named files.
  dframe <- make.theory.dataframe(files)
  if (verbose)
     cat("Made the data frame for distances\n")

  # If we have a "h" column, we want to scale it. We do it here because
  # this is the function that knows about the nature of the directory
  # that it is working on.
  if (any(names(dframe) %in% "h"))
    dframe$h <- dframe$h * 2.99792458e+05

  # Plot each interesting column against "z".
  y.cols <- Filter(function(x) x!="z", names(dframe))
  for (y.col in y.cols)
    make.distance.plot(y.col, dframe, prefix, outdir, verbose, devname)
}

#' Make all the CMB spectrum plots for files under the given top-level
#' directory *topdir*. CMB spectrum plots are made for all data files in
#' the subdirectory "cmb_cl".
make.cmb.spectrum.plots <- function(topdir, verbose, prefix, outdir, devname)
{
  # data files are in "cmb_cl"
  datadir <- file.path(topdir, "cmb_cl")
  if (verbose) cat("Making CMB spectrum plots from data in", datadir, "\n")

  files <- dir(datadir, full.names=TRUE)

  # Create a dataframe from all the named files.
  dframe <- make.theory.dataframe(files)
  if (verbose)
    cat("Made the data frame for CMS spectrum plots\n")

  # Plot each interesting column against "C_ell"
  y.cols <- Filter(function(x) x!="ell", names(dframe))
  for (y.col in y.cols)
    make.cmb.spectrum.plot(y.col, dframe, prefix, outdir, verbose, devname)

  # Make a combined plot of all spectra versus "ell".
  molten<-melt(dframe, id="ell")
  make.cmb.grand.spectrum.plot(molten, prefix, outdir, verbose, devname)
}

make.matter.power.plots <- function(topdir, verbose, prefix, outdir, devname)
{
  # data files are in "matter_power_lin" and "matter_power_nl".
  lin_data_dir <- file.path(topdir, "matter_power_lin")
  nl_data_dir <- file.path(topdir, "matter_power_nl")

  if (verbose)
    cat("Making matter power plots from data in", lin_data_dir, "and", nl_data_dir)

  # We don't have to search what files are present, because the matter
  # power output format is special. Our handling has to rely on that.
  df.lin <- make.matterpower.dataframe(lin_data_dir, "linear")
  df.nl <- make.matterpower.dataframe(nl_data_dir, "nonlinear")

  # Combine the dataframes, creating appropriate factors.
  dframe <- rbind(df.lin, df.nl)
  dframe$type <- as.factor(dframe$type)

  plotfile <- paste(outdir, "/", prefix, "_matter_power.", devname, sep="")
  if (verbose)
    cat("Plotting matter power plot into", plotfile, "\n")
  log.log.scales=list( x=list(log=10, equispaced=FALSE)
                     , y=list(log=10, equispaced=FALSE)
                     )
  device.fcn = get(devname)
  device.fcn(plotfile)
  print( xyplot( p_k~k_h
               , subset(dframe, z==0)
               , ylab="p_k"
               , xlab="k"
               , group=type
               , type="l"
               , grid=TRUE
               , auto.key=list(space="top", lines=TRUE, points=FALSE)
               , scales = scales.log.log()
               , xscale.components = xscale.log()
               , yscale.components = yscale.log()
               ))
  invisible(dev.off())
}

make.data.frame <- function(fname)
{
  d <- read.table(fname)
  first <- readLines(fname, n=1)
  first <- sub("#", "", first)          # Remove comment
  parts <- strsplit(first, "\t")[[1]]   # split on tabs
  cols <- sub("[a-zA-Z_]+--", "", parts) # remove leading section names

  names(d) <- cols
  likes <- exp(d$LIKE)
  norm <- sum(likes)
  d$l <- likes/norm
  return(d)
}

#' Make a 1-d posterior density plot for each variable in the dataframe.
#' We use bw="nrd", which gives a bandwidth calculation according to:
#'    Scott, D. W. (1992) Multivariate Density Estimation: Theory, Practice, and
#'    Visualization. Wiley.
make.1d.likelihood.plots <- function(dframe, prefix, output, device, verbose)
{
  cols <- Filter(function(n) { ! n %in% c("l","LIKE") }, names(dframe))
  for(col in cols)
  {
    if (verbose) cat("Making 1-d density plot for", col, "\n")
    form <- as.formula(paste("l~", col, sep=""))
    dframe.summary <- summaryBy(form, data = dframe, FUN=sum)
    names(dframe.summary)[2] <- "l" # Replace ugly name generated by summaryBy
    dev.fcn <- get(device)
    filename <- file.path(output, paste(prefix, "_", col, ".", device, sep=""))
    p1 <- xyplot( form
                , dframe.summary
                , type="l"
                , lwd=2
                , ylab = "likelihood"
                , grid = TRUE
                )
    dev.fcn(filename)
    print(p1)
    invisible(dev.off())
  }
}

#' Return the likelihood values corresponding to the boundary of the
#' regions containing the probability contents 'levels'.
find.contours <- function(dframe, levels = c(0.68, 0.95))
{
  probs.sorted <- sort(dframe$l, decreasing=TRUE)
  probs.cs <- cumsum(probs.sorted)
  # Get the indices of the first values greater than the given confidence levels.
  indices <- sapply( levels
                   , function(x) which(probs.cs>x)[1]
                   )
  probs.sorted[indices]
}

#' vmat2df will convert the kind of list returned by kde2d (containing two)
#' vectors and a matrix, named x, y and z) into a dataframe with columns
#' x, y, z.
vmat2df <- function(u)
{
  g = expand.grid(u$x, u$y)
  data.frame( x = g$Var1, y = g$Var2, z = as.numeric(u$z))
}

#' Make a 2-d density plot of  xcol vs. ycol, using data from df.
make.2d.density.plot <- function( df, xcol, ycol, prefix, output, device
                                , use.color
                                )
{
  form <- as.formula(paste("l~", xcol, "+", ycol, sep=""))
  df.summary <- summaryBy(form, data = df, FUN = sum)
  names(df.summary)[3] <- "l" # replace the ugly name given by summaryBy

  # Find the values of z which correspond to the given confidence levels.
  conf.levels = c(0.68, 0.95)
  zvals = find.contours(df.summary, conf.levels)
  
  dev.fcn <- get(device)
  filename <- file.path( output
                       , paste(prefix, "_", xcol, "_", ycol, ".", device, sep=""))
  levels <- c(0, zvals, 1)
  labels <- as.character(c(0, conf.levels, 1))
  form <- as.formula(paste("l~", xcol, "*", ycol, sep=""))
  p <- contourplot( form, df.summary, at=levels, labels=labels
                  , panel=function(...){panel.grid(-1,-1); panel.contourplot(...)}
                  , xlab = xcol
                  , ylab = ycol
                  , region = use.color
                  , lwd=2
                  , col.regions = function(n,a) rev(terrain.colors(n,a))
                  , colorkey = FALSE
                  )
  dev.fcn(filename)
  print(p)
  invisible(dev.off())
}

#' For each pair, generate the kde2d result matrix.
#' Determine the values of z at which the 68% and 95% contour lines lie.
#' Transform the matrix to a dataframe.
#' Make the contour plot.
make.2d.density.plots <- function(df, prefix, output, device, verbose, use.color)
{
  if (verbose) cat("Making 2-d density plots\n")
  # Go through all pairs of interesting variables (all but 'LIKE' and 'l',
  # the last two columns).
  n.interesting <- ncol(df)-2
  cols <- names(df)[1:n.interesting]
  pairs <- combn(cols, 2, simplify=FALSE)
  for(pair in pairs)
  {
    xcol <- pair[[1]]
    ycol <- pair[[2]]
    if (verbose) cat("Making 2-d density plot of", xcol, "vs.", ycol, "\n")
    make.2d.density.plot(df, xcol, ycol, prefix, output, device, use.color)
  }

}

#' Create a dataframe from CosmoSIS output. We expect the first line of
#' the output to contain the names of the parameters, separated by
#' spaces, and with section names separated from parameter names by a
#' double-hyphen.
make.data.frame <- function(fname, burn)
{
  d <- read.table(fname)
  # Remove the first 'burn' elements
  d <- d[-c(1:burn),]
  first <- readLines(fname, n=1)
  first <- sub("#", "", first)          # Remove comment
  parts <- strsplit(first, "\t")[[1]]   # split on tabs
  cols <- sub("[a-zA-Z_]+--", "", parts) # remove leading section names
  names(d) <- cols
  likes <- exp(d$LIKE)
  norm <- sum(likes)
  d$l <- likes/norm
  return(d)
}

#' Make a 1-d posterior density plot for each variable in the dataframe.
#' We use bw="nrd", which gives a bandwidth calculation according to:
#'    Scott, D. W. (1992) Multivariate Density Estimation: Theory, Practice, and
#'    Visualization. Wiley.
make.1d.post.density.plots <- function(dframe, prefix, output, device, verbose)
{
  cols <- Filter(function(n) { ! n %in% c("l","LIKE") }, names(dframe))
  for(col in cols)
  {
    if (verbose) cat("Making 1-d density plot for", col, "\n")
    dev.fcn <- get(device)
    filename <- file.path(output, paste(prefix, "_", col, ".", device, sep=""))
    form = as.formula(paste("~", col, sep=""))
    dev.fcn(filename)
    print(densityplot( form
                     , dframe
                     , type="l"
                     , lwd=2
                     , panel = function(...)
                     {
                       panel.grid(-1,-1)
                       panel.densityplot(...)
                     }
                     , ylab = "posterior probability density"
                     , plot.points = FALSE
                     , bw = "nrd"
                     ))
    invisible(dev.off())
  }
}

find.contours <- function(kde, levels = c(0.68, 0.95))
{
  probs.sorted <- sort(kde$z, decreasing=TRUE)
  probs.cs <- cumsum(probs.sorted)
  # Get the indices of the first values greater than the given confidence levels.
  indices <- sapply( levels
                   , function(x) which(probs.cs>x)[1]
                   )
  probs.sorted[indices]
}

#' vmat2df will convert the kind of list returned by kde2d (containing two)
#' vectors and a matrix, named x, y and z) into a dataframe with columns
#' x, y, z.
vmat2df <- function(u)
{
  g = expand.grid(u$x, u$y)
  data.frame( x = g$Var1, y = g$Var2, z = as.numeric(u$z))
}

#' Make a 2-d density plot of  xcol vs. ycol, using data from df.
make.2d.density.plot <- function( df, xcol, ycol, prefix, output, device
                                , use.color
																, nbins)
{
  kde <- kde2d(df[,xcol], df[,ycol], n=nbins	)

  # kde$x and kde$y carry the bin co-ordinates.
  # kde$z carries the estimated density at (x,y). We convert kde$z into
  # the probability for the corresponding bin. We could multiply the density
  # by the area of the bin, but encounter less trouble from rounding if we
  # normalize the sum to 1.
  norm <- sum(kde$z)
  kde$z <- kde$z/norm

  # Find the values of z which correspond to the given confidence levels.
  conf.levels = c(0.68, 0.95)
  zvals = find.contours(kde, conf.levels)
  
  # Convert from vectors+matrix to dataframe, to use lattice plotting.
  d <- vmat2df(kde)
  
  dev.fcn <- get(device)
  filename <- file.path( output
                       , paste(prefix, "_", xcol, "_", ycol, ".", device, sep=""))
  levels <- c(0, zvals, 1)
  labels <- as.character(c(0, conf.levels, 1))
  p <- contourplot( z~x*y, d, at=levels, labels=labels
                  , panel=function(...){panel.grid(-1,-1); panel.contourplot(...)}
                  , xlab = xcol
                  , ylab = ycol
                  , region = use.color
                  , lwd=2
                  , col.regions = function(n,a) rev(terrain.colors(n,a))
                  , colorkey = FALSE
                  )
  dev.fcn(filename)
  print(p)
  invisible(dev.off())
}

#' For each pair, generate the kde2d result matrix.
#' Determine the values of z at which the 68% and 95% contour lines lie.
#' Transform the matrix to a dataframe.
#' Make the contour plot.
make.2d.density.plots <- function(df, prefix, output, device, verbose, use.color, nbins)
{
  # Go through all pairs of interesting variables (all but 'LIKE' and 'l',
  # the last two columns).
  n.interesting <- ncol(df)-2
  cols <- names(df)[1:n.interesting]
  pairs <- combn(cols, 2, simplify=FALSE)
  for(pair in pairs)
  {
    xcol <- pair[[1]]
    ycol <- pair[[2]]
    if (verbose) cat("Making 2-d density plot of", xcol, "vs.", ycol, "\n")
    make.2d.density.plot(df, xcol, ycol, prefix, output, device, use.color, nbins)
  }

}

#' Load data from all the files in fnames. These are expected to be
#' cosmosis timing data files, with the canonical naming pattern:
#'    cosmosis-timing-<pid>.db
#'
#' Returns a data.frame.
#'
load.cosmosis.timing.data <- function(fnames)
{
  drv <- dbDriver("SQLite")
  on.exit(dbUnloadDriver(drv))
  result <- Reduce(function(dframe, fname) {
                    newframe <- load.data.aux(drv,fname)
                    rbind(newframe, dframe)
                    },
                    fnames,
                    data.frame())
  result$pid <- as.factor(result$pid)
  result$module <- as.factor(result$module)
  result
}


load.data.aux <- function(driver, fname)
{
  con <- dbConnect(driver, fname)
  on.exit(dbDisconnect(con))
  tmp <- dbReadTable(con, "Samples")
  tmp$pid <- getPid(fname)
  tmp
}

#' Get the process id from a canonical cosmosis-timing filename. The
#' filename is expected to be in the form:
#'    cosmosis-timing-<pid>.db
#'
getPid <- function(fname)
{
   parts = strsplit(fname, "-|\\.")[[1]]
   if (length(parts) != 3) return(1)
   as.integer(parts[[3]])
}
