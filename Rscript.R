
# --------------------------------------------------------------------------- #
#                               math project final                            #
# --------------------------------------------------------------------------- #


rm(list = ls())
gc()
library(ggplot2) 
library(reshape2) 
library(deSolve) 
library(minpack.lm)

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

# --------------------------------------------------------------------------- #
#                                   rain                                      #
# --------------------------------------------------------------------------- #

rain <- read.csv('./data/pr_1901_2020_ZAF.csv') # load data
names(rain)[1] <- 'rain' # assign name
rain_mean <- tapply(rain$rain, rain$Year, mean)
rain_sd <- tapply(rain$rain, rain$Year, sd)
rain <- data.frame(year = as.numeric(names(rain_mean)), rain = rain_mean, rain_sd = rain_sd)
plot_lowess(rain$year, rain$rain, f = .75, main = 'yearly rain', ax = rain$year)


# --------------------------------------------------------------------------- #
#                                temperature                                  #
# --------------------------------------------------------------------------- #

temp <- read.csv('./data/tas_1901_2020_ZAF.csv')
names(temp)[1] <- 'temp'
temp_mean <-tapply(temp$temp, temp$Year, mean)
temp_sd <-tapply(temp$temp, temp$Year, sd)
temp <- data.frame(year = names(temp_mean), temp = temp_mean, temp_sd = temp_sd)
plot_lowess(temp$year, temp$temp, f = .75, main = 'yearly temp', ax = temp$year)

# --------------------------------------------------------------------------- #
#                               co2 emissions                                 #
# --------------------------------------------------------------------------- #


co2 <- read.csv("./data/co2_raw.txt")
co2 <- co2[co2$country == 'South Africa',c('year','co2')]
co2$co2_change <- c(NA, diff(co2$co2))
co2 <- co2[co2$year > 1901,]
par(mfrow = c(1,2))
plot_lowess(co2$year, co2$co2, f = .75, main = 'yearly co2 cumulative', ax = co2$year)
plot_lowess(co2$year, co2$co2_change, f = .75, main = 'yearly co2 change', ax = co2$year)
par(mfrow = c(1,1))

# --------------------------------------------------------------------------- #
#                            join data together                               #
# --------------------------------------------------------------------------- #

X <- merge(merge(rain, temp, by = 'year'), co2, by = 'year') # join data
X <- X[, c('year','rain','temp','co2')]
X$year <- as.numeric(X$year)
pairs(X)
cor(X)

X <- X[,-4] # remove co2 as highly correlated with temp

# --------------------------------------------------------------------------- #
#                            load count data                                  #
# --------------------------------------------------------------------------- #

dat0 <- rbind(
  read.csv("~/Dropbox/git_ACM40960/data/occurrence1.txt", row.names=1),
  read.csv("~/Dropbox/git_ACM40960/data/occurrence2.txt", row.names=1),
  read.csv("~/Dropbox/git_ACM40960/data/occurrence3.txt", row.names=1),
  read.csv("~/Dropbox/git_ACM40960/data/occurrence4.txt", row.names=1)
)
  
  
dat <- dat0

# remove empty species
dat <- dat[ - which(dat$family == ''), ] 

# sanity check on year
sort(unique(dat$year))
dat <- dat[dat$year > 1900,]
years <- cbind( year = min(dat$year) : max(dat$year) )

# pre process data: create Y such that we will have the individual count
# for each species for each year (min #rows required -> 30)
Y <- data.frame()
min.row <- 30
fams <- unique(dat$family)
for(f in fams){
  i.dat <- dat[dat$family == f, ]
  i.years <- length(unique(i.dat$year))
  if(i.years >= min.row){ # if enough obs
    i.x <- tapply(i.dat$individualCount, i.dat$year, mean) # compute annual mean
    i.x <- cbind(year = names(i.x), individualCount = i.x, species = f) # data frame
    i.x <- merge(years, i.x, by = 'year', all.x = 1) # left join with years
    i.x$species[is.na(i.x$species)] <- f
    Y <- rbind(Y, i.x) # append to main data set 
  }
}

# create Y_wide such that we will have one species per column 
# and the related individual count per row
Y_wide <- tidyr::pivot_wider(Y,
                             names_from = 'species',
                             values_from = 'individualCount'
                             )

# quick visualization
matplot(Y_wide[,1], Y_wide[,-1], type = 'l')

# we observe some "unusual spike" according to a population evolution
# we decide to remove it from the data
to_investigate <- Y[which.max(Y$individualCount), 'species']
matplot(Y_wide[,1], Y_wide[,to_investigate], type = 'l')

# remove Glareolidae
Y <- Y[- which(Y$species == 'Glareolidae'),]
Y_wide <- Y_wide[, - which(colnames(Y_wide) == 'Glareolidae')]

# quick visualization
Y_wide <- data.frame(Y_wide)
matplot(Y_wide[,1], Y_wide[,-1], type = 'l', xaxt = 'n')
axis(1, at = Y_wide[,1], labels = Y_wide[,1])
title('Species count vs time')

# stop('2013 20210722')

# --------------------------------------------------------------------------- #
#                            model of interest                                #
# --------------------------------------------------------------------------- #

# combine data ----------------------------------------------------------------

Y <- merge(Y, X, by = 'year', all.x = TRUE) # join data by year
Y <- Y[order(Y$species), ] # sort data by species
Y$individualCount <- round(
  as.numeric(Y$individualCount)
  ) # cast into num

# inspect individual count distribution: we can try a Pois distr
par(mfrow = c(1,2))
hist(Y$individualCount)
hist(Y$individualCount[Y$individualCount < 50])

# center data -----------------------------------------------------------------
scale <- FALSE
Y$syear <- scale(Y$year, scale = scale)
# Y$stemp <- scale(Y$temp, scale = scale)
# Y$srain <- scale(Y$rain, scale = scale)
# Y$individualCount <- round(as.numeric(Y$individualCount))
# Y <- Y[!is.na(Y$temp),]

# fit model -------------------------------------------------------------------

# define degree
q <- 2

library(gam)
fit_year <- lme4::glmer(
  individualCount ~ ns(syear, q) + (1 | species), 
  family = 'poisson', data = Y
)
summary(fit_year)

# what is the effect of year
t <- sort(unique(Y$syear))
t_pred <- ns(t, q) %*% as.numeric(coef(fit_year)[[1]][1,-1])
plot(t, t_pred, type = 'l', main = 'effect of time')

# inspect residuals: we can see some structure not explained by syear
# therefore we will include further predictors
plot(fitted(fit_year), residuals(fit_year))

# further pred
q_temp <- 2
Y$srain <- scale(Y$rain, scale = TRUE)
Y$stemp <- scale(Y$temp, scale = TRUE)
fit_year_temp <- lme4::glmer(
  # this model returns the best AIC among those explored
  individualCount ~ ns(syear, q_temp) * stemp + (1 | species), 
  family = 'poisson', data = Y
)
summary(fit_year_temp)

# less pronounced but still some structure 
par(mfrow = c(1,2))
plot(fitted(fit_year), residuals(fit_year))
plot(fitted(fit_year_temp), residuals(fit_year_temp))

# stop('2146 20210722')

# species count subset
to_sub <- names(
  sort(table(Y[!is.na(Y$individualCount),'species']), decreasing = TRUE)[1]
)
Y1 <- Y[Y$species == to_sub & complete.cases(Y),] # remove NA obs


# --------------------------------------------------------------------------- #
#                      model rain data with differential eq                   #
# --------------------------------------------------------------------------- #

# rain

rainchange <- function(t, state, parms) {
  
  # browser()

    with(as.list(c(state,parms)), {

      # alpha <- as.numeric(parms["alpha"])
      # beta <- as.numeric(parms["beta"])
      gamma <- as.numeric(parms["gamma"])

      # dR = alpha * rain + beta * temp
      # dR = (alpha + beta * temp) * rain
      dT = gamma * temp

      return(list(c(dT)))
  })
}

stop('0839 20210723')

t <- Y1$time <- Y1$year - min(Y1$year)
yini <- c(rain = Y1$rain[1], temp = Y1$temp[1])
yini <- c(temp = Y1$temp[1])
parms <- c(alpha = -1e-2, beta = 1e-2, gamma = 1e-3)
out <- deSolve::ode(y = yini, time = t, func = rainchange, parms = as.list(parms))
out
plot(out)

# function that calculates residual sum of squares
ssq_rain <- function(parms){

  # browser()

  # inital concentration
  cinit <- c(rain = Y1$rain[1], temp = Y1$temp[1])
  cinit <- c(temp = Y1$temp[1])
  
  # time points for which conc is reported
  t <- sort( unique( c(seq(0,5,0.1), t ) ))

  # solve ODE for a given set of parameters
  out=ode(y=cinit,times=t,func=rainchange,parms=as.list(parms))

  # Filter data that contains time points where data is available
  outdf=data.frame(out)
  outdf=outdf[outdf$time %in% Y1$time,]
  # Evaluate predicted vs experimental residual
  #ssqres <- c(outdf$rain-Y1$rain, outdf$temp - Y1$temp)
  ssqres <- outdf$temp - Y1$temp

  return(ssqres) # return predicted vs experimental residual
}

# parameter fitting using levenberg marquart algorithm
parms_rain <- c(alpha = -1e-2, beta = 1e-2, gamma = 1e-3)
parms_rain <- c(gamma = 1e-3)

# fitting
fitval_rain <- nls.lm(par = parms_rain, fn = ssq_rain)
summary(fitval_rain)

# estimate uncertainty
se <- summary(fitval_rain)$coefficients[1, 'Std. Error']
parest_rain <- c(fitval_rain$par - qnorm(.975) * se,
                 fitval_rain$par,
                 fitval_rain$par + qnorm(.975) * se)

# estimate scenarios
# cinit <- c(rain = Y1$rain[1], temp = Y1$temp[1])
cinit <- c(temp = Y1$temp[1])
out_rain <- ode(y=cinit,times=t,func=rainchange,parms=as.list(parest_rain[2]))
out_rain <- data.frame(out_rain)
# 
# propose longer data 
t_long <- 0:60
Y1 <- merge(Y1, data.frame(time = t_long), all.y = TRUE)
out_rain_long <- ode(y=cinit,times=t_long,func=rainchange,parms=as.list(parest_rain))
out_rain_long <- data.frame(out_rain_long)

# predict longer data
out_rain_lb = ode(y=cinit,times=t_long,func=rainchange,parms=as.list(parest_rain[1]))
out_rain_pe = ode(y=cinit,times=t_long,func=rainchange,parms=as.list(parest_rain[2]))
out_rain_ub = ode(y=cinit,times=t_long,func=rainchange,parms=as.list(parest_rain[3]))
out_rain_par <- merge(merge(out_rain_lb, out_rain_pe, by = 'time'), out_rain_ub, by = 'time')
names(out_rain_par) <- c('time', 'lb', 'pe', 'ub')

out_rain_par_only <- out_rain_par[out_rain_par$time > max(out_rain$time),]

par(mfrow = c(1,2))
# plot(Y1$time, Y1$rain, type = 'b')
# points(out_rain$time, out_rain$rain, col = 'red')

ylim <- range(out_rain_par[,-1])
plot(Y1$time, Y1$temp, type = 'l', ylim = ylim, xaxt = 'n')
lines(out_rain$time, out_rain$temp, col = 'red')
axis(1, at = Y1$time, labels = Y1$year)

plot(t_long, Y1$temp, type = 'l', ylim = ylim, xaxt = 'n')
lines(out_rain_par$time, out_rain_par$pe, col = 'red')
polygon(c(out_rain_par_only$time, rev(out_rain_par_only$time)),
        c(out_rain_par_only$lb, rev(out_rain_par_only$ub)),
        col = ggplot2::alpha('red', .2), border = NA)
axis(1, at = Y1$time, labels = Y1$time + 1975)

# --------------------------------------------------------------------------- #
#                   propose a model for rain evolution by temp                #
# --------------------------------------------------------------------------- #

# need to understand how temp evolves
rain_temp <- merge(rain, temp, by = 'year')
fit_rain <- lm(rain ~ temp, data = rain_temp)

temp_seq <- seq(min(rain_temp$temp), max(rain_temp$temp), length.out = 100)
fitteds <- predict(fit_rain, newdata = data.frame(temp = temp_seq))

par(mfrow = c(1,1))
plot(rain_temp$temp, rain_temp$rain)
lines(temp_seq, fitteds, lwd = 2, col = 'red')




# 
# 
# # --------------------------------------------------------------------------- #
# #                      model temp data with differential eq                   #
# # --------------------------------------------------------------------------- #
# 
# 
# # rain
# tempchange <- function(t,c,parms){# rate constant passed through a list called parms
#   
#   beta <- parms$beta
#   
#   # derivatives dr/dt are computed below
#   r <- rep(0,length(c))
#   r[1] <- beta * c["temp"]
#   
#   return(list(r))
# }
# 
# # function that calculates residual sum of squares
# ssq_temp <- function(parms){
#   #browser()
#   
#   # inital concentration
#   cinit <- c(temp = Y1$temp[1])
#   t=c(seq(0,5,1),Y1$time)
#   t=sort(unique(t))
#   
#   # solve ODE for a given set of parameters
#   out=ode(y=cinit,times=t,func=tempchange,parms=as.list(parms))
#   
#   # Filter data that contains time points where data is available
#   outdf=data.frame(out)
#   outdf=outdf[outdf$time %in% Y1$time,]
#   ssqres=outdf$temp-Y1$temp
#   
#   return(ssqres) # return predicted vs experimental residual
# }
# 
# # parameter fitting using levenberg marquart algorithm
# cinit <- c(temp = Y1$temp[1])
# t <- Y1$time
# parms_temp <- c(beta = 0.1)
# 
# # fitting
# fitval_temp <- nls.lm(par = parms_temp, fn = ssq_temp)
# summary(fitval_temp)
# 
# parest_temp <- fitval_temp$par # parest
# out_temp = ode(y=cinit,times=t_long,func=tempchange,parms=as.list(parest_temp))
# 
# # visualize results ------------------------------------------------------------
# 
# par(las = 2, mfrow = c(1,2))
# 
# plot(out_rain[,1], out_rain[,2], xaxt = 'n', lwd = 2, col = 'royalblue', type = 'l',
#      xlab = 'time', ylab = 'rain')
# title('Evolution of rain in time')
# points(Y1$time, Y1$rain)
# axis(1, at = t_long, t_long + 1975)
# 
# 
# plot(out_temp[,1], out_temp[,2], xaxt = 'n', lwd = 2, col = 'darkorange2', type = 'l',
#      xlab = 'time', ylab = 'temperature')
# title('Evolution of temperature in time')
# points(Y1$time, Y1$temp)
# axis(1, at = t_long, t_long + 1975)
# 
# 
# 
# # i want to estimate the uncertainty of these estimation
#  
# warning('do you really want to specify your matrix in this way?')
# S <- diag(c(vcov(fitval_rain), vcov(fitval_temp)))
# dof <- nrow(Y1) * 3 #dof
# 
# # draw the confidence region
# # get points for a circle with radius r
# r=sqrt(qf(0.95,2,dof)*2)
# theta=seq(0,2*pi,length.out=100)
# z=cbind(r*cos(theta),r*sin(theta))
# # transform points of circle into points of ellipse using
# # svd of covariance matrix
# S_svd=svd(S)      # SVD of covariance matrix
# xt=t(S_svd$v)%*%diag(sqrt(S_svd$d))%*%t(z) # transform from circle to ellispse
# x=t(xt)
# # translate the ellipse so that center is the estimated parameter value
# parest <- c(alpha = parest_rain, beta = parest_temp)
# x=x+matrix(rep(as.numeric(parest),100),nrow=100,byrow=T)
# 
# par(mfrow = c(1,1))
# plot(x[,1],x[,2],type="l",xlab="alpha",ylab="beta",lwd=2)
# points(parest[1], parest[2], pch=20,col="blue",cex=2)



















