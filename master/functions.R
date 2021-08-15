# functions.R

# --------------------------------------------------------------------------- #
#                               useful functions                              #
# --------------------------------------------------------------------------- #

plot_lowess <- function(
  x, # input data
  y, # response
  f = .5, # span value
  fit.col = 'violetred', # fit line color
  ax = NULL,  # x axis values to be printed on the plot 
  add = 0, # 0: make new plot, 1: superimpose on existing last plot 
  ... # other optional parameter for the plot() function
  ){
  
  # fit non parametric model gien a span value
  model <- lowess(x, y, f = f) # estimate local weighted regression
  fits <- model$y # extract fitted values

  # produce plot of the data and fitted model
  par(las = 2)
  if(!add) plot(x, y, xaxt = 'n', col = 'grey70', ...) # if: create new plot
  if(!is.null(ax)) axis(1, at = x, labels = ax, tick = FALSE) # if: add axis values
  lines(x, fits, lwd = 2, col = fit.col) # superimpose fitted values line
  
  # return object only if assigned 
  invisible(fits)
}
