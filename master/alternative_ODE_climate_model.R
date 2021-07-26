
# --------------------------------------------------------------------------- #
#                               math project final                            #
# --------------------------------------------------------------------------- #


rm(list = ls()) # clean env
source('functions.R') # load functions  

gc()
library(ggplot2) 
library(reshape2) 
library(deSolve) 
library(minpack.lm)


# --------------------------------------------------------------------------- #
#                                   rain                                      #
# --------------------------------------------------------------------------- #

rain <- read.csv('../data/pr_1901_2020_ZAF.csv') # load data
names(rain)[1] <- 'rain' # assign name
rain_mean <- tapply(rain$rain, rain$Year, mean)
rain_sd <- tapply(rain$rain, rain$Year, sd)
rain <- data.frame(year = as.numeric(names(rain_mean)), rain = rain_mean, rain_sd = rain_sd)
plot_lowess(rain$year, rain$rain, f = .75, main = 'yearly rain', ax = rain$year)


# --------------------------------------------------------------------------- #
#                                temperature                                  #
# --------------------------------------------------------------------------- #

temp <- read.csv('../data/tas_1901_2020_ZAF.csv')
names(temp)[1] <- 'temp'
temp_mean <-tapply(temp$temp, temp$Year, mean)
temp_sd <-tapply(temp$temp, temp$Year, sd)
temp <- data.frame(year = names(temp_mean), temp = temp_mean, temp_sd = temp_sd)
plot_lowess(temp$year, temp$temp, f = .75, main = 'yearly temp', ax = temp$year)

# --------------------------------------------------------------------------- #
#                               co2 emissions                                 #
# --------------------------------------------------------------------------- #


# co2 <- read.csv("./data/co2_raw.txt")
# co2 <- co2[co2$country == 'South Africa',c('year','co2')]
# co2$co2_change <- c(NA, diff(co2$co2))
# co2 <- co2[co2$year > 1901,]
# par(mfrow = c(1,2))
# plot_lowess(co2$year, co2$co2, f = .75, main = 'yearly co2 cumulative', ax = co2$year)
# plot_lowess(co2$year, co2$co2_change, f = .75, main = 'yearly co2 change', ax = co2$year)
# par(mfrow = c(1,1))

# --------------------------------------------------------------------------- #
#                            join data together                               #
# --------------------------------------------------------------------------- #


# X <- merge(merge(rain, temp, by = 'year'), co2, by = 'year') # join data
X <- merge(rain, temp, by = 'year')
# X <- X[, c('year','rain','temp','co2')]
X <- X[, c('year','rain','temp')]
pairs(X)
cor(X)

# X <- X[,-4] # remove co2 as highly correlated with temp


# --------------------------------------------------------------------------- #
#                      model rain data with differential eq                   #
# --------------------------------------------------------------------------- #

warning('here it would be intersting to track the evolution of temp and rain together')

# rain

rainchange <- function(t, state, parms) {
  
  # browser()
  
  with(as.list(c(state,parms)), {
    
    #alpha <- as.numeric(parms["alpha"])
    beta <- as.numeric(parms["beta"])
    delta <- as.numeric(parms["delta"])
    gamma <- as.numeric(parms["gamma"])
    
    dR = (beta * temp + delta)*rain
    dT = gamma * temp 
#    dT = gamma * temp + alpha * rain
    
    return(list(c(dR, dT)))
  })
}

t <- X$time <- X$year - min(X$year, na.rm = TRUE)
yini <- c(rain = mean(X$rain[1:5]), temp = mean(X$temp[1:5]))
parms <- c(beta = -1e-2, gamma = 1e-3, delta = -1e-4) #, alpha = 1e-4)
out <- deSolve::ode(y = yini, time = t, func = rainchange, parms = as.list(parms))
plot(out)

# function that calculates residual sum of squares
ssq_rain <- function(parms){
  
  #browser()
  
  # inital concentration
  cinit <- c(rain = mean(X$rain[1:5]), temp = mean(X$temp[1:5]))
  
  # time points for which conc is reported
  t <- sort( unique( c(seq(0,5,0.1), t ) ))
  
  # solve ODE for a given set of parameters
  out=ode(y=cinit,times=t,func=rainchange,parms=as.list(parms))
  
  # Filter data that contains time points where data is available
  outdf=data.frame(out)
  outdf=outdf[outdf$time %in% X$time,]
  # Evaluate predicted vs experimental residual
  #ssqres <- c(outdf$rain-Y1$rain, outdf$temp - Y1$temp)
  ssqres <- c(outdf$temp - X$temp,
              outdf$rain - X$rain
  )
  
  return(ssqres) # return predicted vs experimental residual
}

# parameter fitting using levenberg marquart algorithm
# parms_rain <- c(beta = -1e-2, gamma = 1e-3, alpha = 1e-4, delta = -1e-4)
parms_rain <- parms

# fitting
rm(fitval_rain)
fitval_rain <- nls.lm(par = parms_rain, fn = ssq_rain, control = list(maxiter = 1e2))
summary(fitval_rain)

# viz result
parms <- fitval_rain$par
outfit <- deSolve::ode(y = yini, time = t, func = rainchange, parms = as.list(parms))

par(mfrow = c(1,2))
ylim <- range(X$rain)
plot(outfit[,1], outfit[,2], type = 'l', ylim = ylim)
points(X$time, X$rain)
ylim <- range(X$temp)
plot(outfit[,1], outfit[,3], type = 'l', ylim = ylim)
points(X$time, X$temp)


# estimate uncertainty
se <- summary(fitval_rain)$coefficients[1, 'Std. Error']
parest_rain <- c(fitval_rain$par - qnorm(.975) * se,
                 fitval_rain$par,
                 fitval_rain$par + qnorm(.975) * se)

# estimate scenarios
# cinit <- c(rain = Y1$rain[1], temp = Y1$temp[1])
cinit <- c(temp = mean(X$temp[1:5]))
out_rain <- ode(y=cinit,times=t,func=rainchange,parms=as.list(parest_rain[2]))
out_rain <- data.frame(out_rain)
# 
# propose longer data 
t_long <- 0:150
X <- merge(X, data.frame(time = t_long), all.y = TRUE)
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

ylim <- range(out_rain_par[,-1], X$temp, na.rm = TRUE)
plot(X$time, X$temp, type = 'l', ylim = ylim, xaxt = 'n')
lines(out_rain$time, out_rain$temp, col = 'red')
axis(1, at = X$time, labels = X$year)

plot(t_long, X$temp, type = 'l', ylim = ylim, xaxt = 'n')
lines(out_rain_par$time, out_rain_par$pe, col = 'red')
polygon(c(out_rain_par_only$time, rev(out_rain_par_only$time)),
        c(out_rain_par_only$lb, rev(out_rain_par_only$ub)),
        col = ggplot2::alpha('red', .2), border = NA)
axis(1, at = X$time, labels = X$time + min(X$year, na.rm = 1))

# stop('after temp model')

# --------------------------------------------------------------------------- #
#                   propose a model for rain evolution by temp                #
# --------------------------------------------------------------------------- #

# need to understand how temp evolves
# rain_temp <- merge(rain, temp, by = 'year')
rain_temp <- X[complete.cases(X),]

hist(rain_temp$rain) # approximately normal

fit_rain <- lm(rain ~ temp, data = rain_temp) # fit model

temp_seq <- seq(min(rain_temp$temp), max(rain_temp$temp), length.out = 100)
fitteds <- predict(fit_rain, newdata = data.frame(temp = temp_seq))

par(mfrow = c(1,1))
plot(rain_temp$temp, rain_temp$rain)
lines(temp_seq, fitteds, lwd = 2, col = 'red')

# model diagnostic: may use weighted least square
par(mfrow = c(1,2))
plot(fitted(fit_rain), residuals(fit_rain), main = 'ordinary')

# extract wi
wi <- 1 / residuals(fit_rain)**2
fit_rain_wi <- lm(rain ~ temp, data = rain_temp, weights = wi)
plot(fitted(fit_rain_wi), residuals(fit_rain_wi), main = 'weighted')
summary(fit_rain_wi)

# predict rain from temp
pred_type <- 'prediction'
rain_pred_pe <- predict(fit_rain_wi, newdata = data.frame(temp = out_rain_par$pe), interval = pred_type)[,'fit']
rain_pred_lb <- predict(fit_rain_wi, newdata = data.frame(temp = out_rain_par$lb), interval = pred_type)[,'lwr']
rain_pred_ub <- predict(fit_rain_wi, newdata = data.frame(temp = out_rain_par$ub), interval = pred_type)[,'upr']

par(mfrow = c(1,2))
# plot temp as well
plot(rain_temp$time, rain_temp$temp, type = 'l', xlim = range(out_rain_par$time), xaxt = 'n')
title('Temp vs Time')
axis(1, at = out_rain_par$time, labels = out_rain_par$time + min(X$year, na.rm = TRUE))
lines(out_rain_par$time, out_rain_par$pe, lwd = 2, col = 'red')
idx <- out_rain_par$time > max(rain_temp$time)
polygon(c(out_rain_par$time[idx], rev(out_rain_par$time[idx])),
        c(out_rain_par$lb[idx], rev(out_rain_par$ub[idx])),
        border = NA, col = ggplot2::alpha('red', .2))



plot(rain_temp$time, rain_temp$rain, type = 'l', xlim = range(out_rain_par$time), xaxt = 'n')
axis(1, at = out_rain_par$time, labels = out_rain_par$time + min(X$year, na.rm = TRUE))
title('Rain vs Time')
lines(out_rain_par$time, rain_pred_pe, lwd = 2, col = 'red')
idx <- out_rain_par$time > max(rain_temp$time)
polygon(c(out_rain_par$time[idx], rev(out_rain_par$time[idx])),
        c(rain_pred_lb[idx], rev(rain_pred_ub[idx])),
        border = NA, col = ggplot2::alpha('red', .2))


# attach new data to X
X_pred <- cbind(X, rain_pred_lb, rain_pred_pe, rain_pred_ub, out_rain_par[,-1])
names(X_pred) <- c('time', 'year', 'rain_obs', 'temp_obs',
                   'rain_pred_lb', 'rain_pred_pe', 'rain_pred_ub',
                   'temp_pred_lb', 'temp_pred_pe', 'temp_pred_ub'
                   )


# write data
# write.csv(X_pred, '../data/climate_pred.csv')

