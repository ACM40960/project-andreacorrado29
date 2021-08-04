
# list all packages
pkgs <- installed.packages()[,1]

# packages for the code
useful <- c('ggplot2', 'reshape2', 'deSolve', 'minpack.lm','gam', 'xtable','e1071')

# install those not available
for (p in useful) if (!(p %in% pkgs)) install.packages(p)
