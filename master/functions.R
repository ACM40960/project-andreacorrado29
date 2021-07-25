# functions.R

# --------------------------------------------------------------------------- #
#                               useful functions                              #
# --------------------------------------------------------------------------- #

plot_lowess <- function(x, y, f = .5, title = '', ax = NULL, type = 'l', ...){
  
  # fit model
  fits <- lowess(x, y, f = f)$y
  ratio <- fits[NROW(fits)]/fits[1]
  
  # produce plot
  par(las = 2)
  plot(x, y, type = type, xaxt = 'n', col = 'grey50', ...)
  title(title)
  if(!is.null(ax)) axis(1, at = x, labels = ax)
  lines(x, fits, lwd = 2, col = 'violetred')
  
  # return object
  invisible(fits)
}
