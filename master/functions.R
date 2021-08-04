# functions.R

# --------------------------------------------------------------------------- #
#                               useful functions                              #
# --------------------------------------------------------------------------- #

plot_lowess <- function(x, y, f = .5, fit.col = 'violetred', ax = NULL, add = 0, ...){
  
  # fit non parametric model gien a span value
  model <- lowess(x, y, f = f) 
  fits <- model$y

  # produce plot of the data and fitted model
  par(las = 2)
  if(!add) plot(x, y, xaxt = 'n', col = 'grey70', ...)
  if(!is.null(ax)) axis(1, at = x, labels = ax, tick = FALSE)
  lines(x, fits, lwd = 2, col = fit.col)
  
  # return object
  invisible(fits)
}
