
# --------------------------------------------------------------------------- #
#                               math project final                            #
# --------------------------------------------------------------------------- #


rm(list = ls()) # clean env
gc() # clean memory
source('functions.R') # load useful functions  


source('pacakges_install.R') # install libreries if not available

# load useful libraries
library(e1071)
library(xtable)
library(ggplot2) 
library(reshape2) 
library(deSolve) 
library(minpack.lm)

# --------------------------------------------------------------------------- #
#                                   rain                                      #
# --------------------------------------------------------------------------- #

rain <- read.csv('../data/pr_1901_2020_ZAF.csv') # load rain data
names(rain)[1] <- 'rain' # assign name to data set
rain_mean <- tapply(rain$rain, rain$Year, mean) # compute annual mean
# collect data in a easy to handle data set
rain <- data.frame(year = as.numeric(names(rain_mean)), rain = rain_mean)
# plot non parametric local regression to the data of interest + graphics parameter

# span value: 1
plot_lowess(rain$year, rain$rain, f = 1, 
            main = 'Average rain fall amount per year', ax = rain$year,
            xlab = 'year', ylab = 'rain fall amount', type = 'l')

# span value: .5
plot_lowess(rain$year, rain$rain, f = .5, 
            main = 'Average rain fall amount per year', ax = rain$year,
            xlab = 'year', ylab = 'rain fall amount', type = 'l',
            add = 1, fit.col = 'royalblue')

# add legend
legend('topright', legend = c('span 100%', 'span  50%'), lty = 1, lwd = 2,
       bty = 'n', col = c('violetred', 'royalblue'))

# --------------------------------------------------------------------------- #
#                                temperature                                  #
# --------------------------------------------------------------------------- #

temp <- read.csv('../data/tas_1901_2020_ZAF.csv') # load temp data
names(temp)[1] <- 'temp' # assign name to data set
temp_mean <-tapply(temp$temp, temp$Year, mean) # compute average year temp
# assign to an handy data set
temp <- data.frame(year = names(temp_mean), temp = temp_mean)
# fit local regression and produce plot of the result

# span 100%
plot_lowess(temp$year, temp$temp, f = 1, ax = temp$year,
            xlab = 'year', ylab = 'temperature', type = 'l',
            main = 'Annual temperature over year')

# span 50%
plot_lowess(temp$year, temp$temp, f = .5, main = 'yearly temp', ax = temp$year,
            xlab = 'year', ylab = 'temperature', type = 'l',
            add = 1, fit.col = 'blue')

# add legend
legend('topleft', legend = c('span 100%', 'span  50%'), lty = 1, lwd = 2,
       bty = 'n', col = c('violetred', 'royalblue'))



# --------------------------------------------------------------------------- #
#                            join data together                               #
# --------------------------------------------------------------------------- #

X <- merge(rain, temp, by = 'year') # join data in a dataset by year 
pairs(X) # shows pairs plot to give an idea of the relationship
cor(X) # output correlation table
xtable(cor(X), type = "latex") # output table for latex

# the dependence of time in over time is very clear -> propose an ODE model
# the dependence of rain over time is somewhat weak, however, we can observe
# a clear structure in the relation btw rain amount and temperature -> lm

# --------------------------------------------------------------------------- #
#                      model rain data with differential eq                   #
# --------------------------------------------------------------------------- #

# model temperature in time with gamma parameter
tempchange <- function(t, state, parms) {
  
  with(as.list(c(state,parms)), {
    
    gamma <- as.numeric(parms["gamma"]) # gamma parameter
    dtemp <- gamma * temp  # temperature evolution dxdt = gamma x 
    
    return(list(c(dtemp))) # return object as a list
  })
}

# numerically integrate ODE
t <- X$time <- X$year - min(X$year, na.rm = TRUE) # extract time from t = 0 for ODE
yini <- c(temp = mean(X$temp[1:5])) # assign initial condition as the average of the first 5y
parms <- c(gamma = 1e-3) # give initial guess for the parametrs

# function that calculates residual sum of squares for LM algorithm
ssq_temp <- function(parms){
  
  # inital temperature as 5y average
  cinit <- c(temp = mean(X$temp[1:5]))
  
  # time points for which temperature is reported
  t <- sort( unique( c(seq(0,5,0.1), t ) ))
  
  # solve ODE for a given set of parameters
  out=ode(y=cinit,times=t,func=tempchange,parms=as.list(parms))
  
  # Filter data that contains time points where data is available
  outdf=data.frame(out)
  outdf=outdf[outdf$time %in% X$time,]
  # Evaluate predicted vs experimental residual
  ssqres <- outdf$temp - X$temp
  
  return(ssqres) # return predicted vs experimental residual
}


# fitting LM algorithm to the residual sum of square returned by the function ssq_temp
fitval_temp <- nls.lm(par = parms, fn = ssq_temp)
summary(fitval_temp) # provide point estimate + standard error

# estimate uncertainty with 95% confidence interval for the parameter of interest
se <- summary(fitval_temp)$coefficients[1, 'Std. Error']
parest_temp <- c(fitval_temp$par - qnorm(.975) * se,
                 fitval_temp$par,
                 fitval_temp$par + qnorm(.975) * se)

# estimate scenarios with point estimate
out_temp <- ode(y=yini,times=t,func=tempchange,parms=as.list(parest_temp[2]))
out_temp_l <- ode(y=yini,times=t,func=tempchange,parms=as.list(parest_temp[1]))
out_temp_u <- ode(y=yini,times=t,func=tempchange,parms=as.list(parest_temp[3]))

# plot point estimate
ylim <- range(X$temp)
plot(out_temp, ylim = ylim, main = 'Annual temperature ODE evolution',
     xlab = 'year', ylab = 'temperature', lwd = 2, col = 'violetred')
points(X$year - min(X$year), X$temp, col = 'grey70')
polygon(c(out_temp_l[,1], rev(out_temp_u[,1])),
        c(out_temp_l[,2], rev(out_temp_u[,2])), 
        border = NA, col = ggplot2::alpha('red', .1)
)
legend('topleft', col = 'violetred', lwd = 2, lty = 1, 
       legend = 'ODE estimation', bty = 'n')


stop('arrived here')

# propose longer data for future prediction
t_long <- 0:170
X <- merge(X, data.frame(time = t_long), all.y = TRUE) # join data with future years
out_temp_long <- ode(y=yini,times=t_long,func=tempchange,parms=as.list(parest_temp))
out_temp_long <- data.frame(out_temp_long)

# predict longer data with lower, upper bound + point estimate
out_temp_lb = ode(y=yini,times=t_long,func=tempchange,parms=as.list(parest_temp[1]))
out_temp_pe = ode(y=yini,times=t_long,func=tempchange,parms=as.list(parest_temp[2]))
out_temp_ub = ode(y=yini,times=t_long,func=tempchange,parms=as.list(parest_temp[3]))
out_temp_par <- merge(merge(out_temp_lb, out_temp_pe, by = 'time'), out_temp_ub, by = 'time')
names(out_temp_par) <- c('time', 'lb', 'pe', 'ub') # assign name to data set 

# subset only out of sample year, those to be predicted
out_temp_par_only <- out_temp_par[out_temp_par$time > max(out_temp[,'time']),]

# plot result
par(mfrow = c(1,1))
ylim <- range(out_temp_par[,-1], X$temp, na.rm = TRUE)
plot(t_long, X$temp, type = 'l', ylim = ylim, xaxt = 'n')
lines(out_temp_par$time, out_temp_par$pe, col = 'red')
polygon(c(out_temp_par_only$time, rev(out_temp_par_only$time)),
        c(out_temp_par_only$lb, rev(out_temp_par_only$ub)),
        col = ggplot2::alpha('red', .2), border = NA)
axis(1, at = X$time, labels = X$time + min(X$year, na.rm = 1))

# stop('after temp model')


# --------------------------------------------------------------------------- #
#               local sensitivity analysis of temperature evolution           #
# --------------------------------------------------------------------------- #


tempsens <- function(t, state, parms) {
  
  with(as.list(c(state,parms)), {
    
    gamma <- as.numeric(parms["gamma"])
    RHS1 <- gamma * temp 
    RHS2 <- gamma * s + temp # derive sensitivity equation
    
    return(list(c(RHS2, RHS2)))
  })
}


# sensitivity: being only one parameter LSA is OK
# very solid estimation up to 100 SE

# increase qs if you want to investigate more extreme scenario

qs <- 10 ** (1:1) 
par(mfrow = c(length(qs), 2))
for(q in qs){
  
  sens <- c(yini, s = 0)
  pe0 <- parest_temp[2]
  pe1 <- pe0 + q * se
  pe2 <- pe0 - q * se
  out_temp_sens1 <- ode(y=sens,times=t,func=tempsens,parms=list(gamma = pe0))
  out_temp_sens2 <- ode(y=sens,times=t,func=tempsens,parms=list(gamma = pe1))
  out_temp_sens3 <- ode(y=sens,times=t,func=tempsens,parms=list(gamma = pe2))
  
  f <- function(x) log(x)
  ylim <- f(range(out_temp_sens1[,2], out_temp_sens2[,2], out_temp_sens3[,2]))
  plot(out_temp_sens1[,1],  f(out_temp_sens1[,2]), type = 'l', xlab = 'time',
       ylab = '', main = paste('temp vs time \\pm', q,'se'), ylim = ylim)
  lines(out_temp_sens2[,1], f(out_temp_sens2[,2]), type = 'l', col = 'red')
  lines(out_temp_sens3[,1], f(out_temp_sens3[,2]), type = 'l', col = 'blue')
  
  ylim <- f(range(out_temp_sens1[,3], out_temp_sens2[,3], out_temp_sens3[,3]) + 1)
  plot(out_temp_sens1[,1],  f(out_temp_sens1[,3]),
       main = paste('wrt gamma \\pm', q, 'se'), ylab = '', 
       xlab = 'time', type = 'l', ylim = ylim)
  lines(out_temp_sens2[,1], f(out_temp_sens2[,3]), col  = 'red')
  lines(out_temp_sens3[,1], f(out_temp_sens3[,3]), col  = 'blue')
}
par(mfrow = c(1,1))


# --------------------------------------------------------------------------- #
#                   propose a model for rain evolution by temp                #
# --------------------------------------------------------------------------- #

# need to understand how temp evolves
# rain_temp <- merge(rain, temp, by = 'year')
rain_temp <- X[complete.cases(X),]

hist(rain_temp$rain) # approximately normal

fit_rain <- lm(rain ~ splines::ns(temp, 2), data = rain_temp) # fit model
BIC(fit_rain)
summary(fit_rain)
par(mfrow = c(1,2))
plot(fitted(fit_rain), residuals(fit_rain), main = 'ordinary')

temp_seq <- seq(min(rain_temp$temp), max(rain_temp$temp), length.out = 100)
fitteds <- predict(fit_rain, newdata = data.frame(temp = temp_seq))

par(mfrow = c(1,2))
plot(rain_temp$time, rain_temp$rain)
plot(rain_temp$temp, rain_temp$rain)
lines(temp_seq, fitteds, lwd = 2, col = 'red')

# model diagnostic: may use weighted least square
par(mfrow = c(1,2))
plot(fitted(fit_rain), residuals(fit_rain), main = 'ordinary')

# extract wi
wi <- 1 / residuals(fit_rain)**2
fit_rain_wi <- lm(rain ~ splines::ns(temp,2), data = rain_temp, weights = wi)
plot(fitted(fit_rain_wi), residuals(fit_rain_wi), main = 'weighted')
summary(fit_rain_wi)

# predict rain from temp
pred_type <- 'prediction'
rain_pred_pe <- predict(fit_rain_wi, newdata = data.frame(temp = out_temp_par$pe), interval = pred_type)[,'fit']
rain_pred_lb <- predict(fit_rain_wi, newdata = data.frame(temp = out_temp_par$lb), interval = pred_type)[,'lwr']
rain_pred_ub <- predict(fit_rain_wi, newdata = data.frame(temp = out_temp_par$ub), interval = pred_type)[,'upr']


# sensitivity
# try to impute different values
q <- qnorm(.995)
summary(fit_rain_wi)
fit_rain_wi_edit <- fit_rain_wi
fit_rain_wi_edit$coefficients[2] <-  -8.962673 + q *  1.0267
fit_rain_wi_edit$coefficients[3] <-  -8.6207 + q *  0.3338
pred_type <- 'prediction'
rain_pred_pe_e <- predict(fit_rain_wi_edit, newdata = data.frame(temp = out_temp_par$pe), interval = pred_type)[,'fit']
rain_pred_lb_e <- predict(fit_rain_wi_edit, newdata = data.frame(temp = out_temp_par$lb), interval = pred_type)[,'lwr']
rain_pred_ub_e <- predict(fit_rain_wi_edit, newdata = data.frame(temp = out_temp_par$ub), interval = pred_type)[,'upr']

fit_rain_wi_edit$coefficients[2] <-  -8.962673 - q *  1.0267
fit_rain_wi_edit$coefficients[3] <-  -8.6207 - q *  0.3338
pred_type <- 'prediction'
rain_pred_pe_f <- predict(fit_rain_wi_edit, newdata = data.frame(temp = out_temp_par$pe), interval = pred_type)[,'fit']
rain_pred_lb_f <- predict(fit_rain_wi_edit, newdata = data.frame(temp = out_temp_par$lb), interval = pred_type)[,'lwr']
rain_pred_ub_f <- predict(fit_rain_wi_edit, newdata = data.frame(temp = out_temp_par$ub), interval = pred_type)[,'upr']

# viz result
par(mfrow = c(1,2))
# plot temp as well
plot(rain_temp$time, rain_temp$temp, type = 'l',
     xlim = range(out_temp_par$time), xaxt = 'n', col = 'grey70')
title('Temp vs Time')
axis(1, at = out_temp_par$time, labels = out_temp_par$time + min(X$year, na.rm = TRUE))
lines(out_temp_par$time, out_temp_par$pe, lwd = 2, col = 'red')
idx <- out_temp_par$time > max(rain_temp$time)
polygon(c(out_temp_par$time[idx], rev(out_temp_par$time[idx])),
        c(out_temp_par$lb[idx], rev(out_temp_par$ub[idx])),
        border = NA, col = ggplot2::alpha('red', .2))
legend('topleft', col = 'red', lwd = 2, bty = 'n',
       legend = 'Estimated par')

plot(rain_temp$time, rain_temp$rain, type = 'l', xlim = range(out_temp_par$time),
     xaxt = 'n', col = 'grey70')
axis(1, at = out_temp_par$time, labels = out_temp_par$time + min(X$year, na.rm = TRUE))
title('Rain vs Time')
lines(out_temp_par$time, rain_pred_pe, lwd = 2, col = 'red')
idx <- out_temp_par$time > max(rain_temp$time)
polygon(c(out_temp_par$time[idx], rev(out_temp_par$time[idx])),
        c(rain_pred_lb[idx], rev(rain_pred_ub[idx])),
        border = NA, col = ggplot2::alpha('red', .2))
# other scenario + se
lines(out_temp_par$time, rain_pred_pe_e, lwd = 2, col = 'blue')
polygon(c(out_temp_par$time[idx], rev(out_temp_par$time[idx])),
        c(rain_pred_lb_e[idx], rev(rain_pred_ub_e[idx])),
        border = NA, col = ggplot2::alpha('blue', .1))
# other scenario - se
lines(out_temp_par$time, rain_pred_pe_f, lwd = 2, col = 'green')
polygon(c(out_temp_par$time[idx], rev(out_temp_par$time[idx])),
        c(rain_pred_lb_f[idx], rev(rain_pred_ub_f[idx])),
        border = NA, col = ggplot2::alpha('green', .1))
legend('topright', col = c('red','blue','green'), bty = 'n', lty = 1, lwd = 2,
       legend = c('Estimated pars',
                  paste('est +', round(q, 2), 'SE'),
                  paste('est -', round(q, 2), 'SE')
                  )
)

# attach new data to X
X_pred <- cbind(X, rain_pred_lb, rain_pred_pe, rain_pred_ub, out_temp_par[,-1])
names(X_pred) <- c('time', 'year', 'rain_obs', 'temp_obs',
                   'rain_pred_lb', 'rain_pred_pe', 'rain_pred_ub',
                   'temp_pred_lb', 'temp_pred_pe', 'temp_pred_ub'
                   )


# write data
# write.csv(X_pred, '../data/climate_pred.csv')

