# final_code.R

# ---------------------------------------------------------------------- #
# --------------environment set up and useful packages ----------------- #
# ---------------------------------------------------------------------- #

rm(list = ls()) # clean env
gc() # clean memory

source('functions.R') # load useful functions  
source('pacakges_install.R') # install libraries if not available

# load useful libraries
library(ggplot2) 
library(reshape2) 
library(deSolve) 
library(minpack.lm)
library(gam)
library(reshape2)
library(deSolve)
library(viridis)
library(xtable)
library(viridisLite)
library(splines)


# ---------------------------------------------------------------------- #
# -------------------- exploratory data analysis ----------------------- #
# ---------------------------------------------------------------------- #


# rain data ------------------------------------------------------------ #


rain <- read.csv('../data/pr_1901_2020_ZAF.csv') # load rain data
names(rain)[1] <- 'rain' # assign name to data set
rain_mean <- tapply(rain$rain, rain$Year, mean) # compute annual mean
rain_sd <- tapply(rain$rain, rain$Year, sd) # compute annual standard deviation

# collect data in a easy to handle data frame
rain <- data.frame( #
    year = as.numeric(names(rain_mean)), # extract years
    rain = rain_mean # annual rainfall mean
    )

# plot histogram
hist(rain$rain, xlab = 'rainfall amount', main = 'Annual rain distribution', 
     col = 'grey90', border = 'white', breaks = 10)


# temp data ------------------------------------------------------------ #


temp <- read.csv('../data/tas_1901_2020_ZAF.csv') # load temp data
names(temp)[1] <- 'temp' # assign name to data set
temp_mean <- tapply(temp$temp, temp$Year, mean) # compute average year temp
temp_sd   <- tapply(temp$temp, temp$Year, sd) # compute average year temp
# assign to an handy data frame
temp <- data.frame(
    year = names(temp_mean), # extract years
    temp = temp_mean,  # attach mean
    temp_sd = temp_sd # attach sd
    )

# plot mean histogram
hist(temp$temp, xlab = 'temperature', main = 'Annual temperature distribution', 
     col = 'grey90', border = 'white', breaks = 10)



# check relationship with temperature sd and mean for rain --------------------

# produce scatterplot and add correlation 
plot(rain_mean, rain_sd, xlab = 'rainfall year average', 
     ylab = 'rainfall year standard deviation', col = 'grey60')
legend('topleft', col = 'grey', lty = NA, lwd = NA, bty = 'n', pch = 1,
       legend = paste('sd and mean cor =', round(cor(rain_mean, rain_sd), 2)))


# check relationship with temperature sd and mean for temp --------------------

plot(temp_mean, temp_sd, xlab = 'temperature annual average', 
     ylab ='rainfall year standard deviation', col = 'grey60')
legend('topright', col = 'grey', lty = NA, lwd = NA, bty = 'n', pch = 1,
       legend = paste('sd and mean cor =', round(cor(temp_mean, temp_sd), 2)))

# check relationship with temperature sd and mean for temp increment ----------

plot(c(0, diff(temp_mean)), temp_sd, xlab = 'temperature annual increment', 
     ylab ='rainfall year standard deviation', col = 'grey60')
legend('topleft', col = 'grey', lty = NA, lwd = NA, bty = 'n', pch = 1,
       legend = paste('sd and mean cor =', round(cor(c(0, diff(temp_mean)), temp_sd), 2)))


# plot non parametric local regression to for rainfall ----------------------------

# span value: 1
plot_lowess(rain$year, rain$rain, f = 1, 
            main = 'Average rain fall amount per year', ax = rain$year,
            xlab = 'year', ylab = 'rain fall amount', type = 'l')
legend('topright', legend = c('span 100%'), lty = 1, lwd = 2,
       bty = 'n', col = c('violetred'))

# span value: .5
plot_lowess(rain$year, rain$rain, f = .5, 
            main = 'Average rain fall amount per year', ax = rain$year,
            xlab = 'year', ylab = 'rain fall amount', type = 'l',
            add = 0, fit.col = 'royalblue')

# add legend
legend('topright', legend = c('span  50%'), lty = 1, lwd = 2,
       bty = 'n', col = c('royalblue'))


# plot non parametric local regression to for temperature ----------------------------


# span 1
plot_lowess(temp$year, temp$temp, f = 1, ax = temp$year,
            xlab = 'year', ylab = 'temperature', type = 'l',
            main = 'Annual temperature over year')

legend('topleft', legend = c('span  100%'), lty = 1, lwd = 2,
       bty = 'n', col = c('violetred'))

# span .5
plot_lowess(temp$year, temp$temp, f = .5, ax = temp$year,
            xlab = 'year', ylab = 'temperature', type = 'l',
            add = 0, fit.col = 'blue', main = 'Annual temperature over year')
legend('topleft', legend = c('span  50%'), lty = 1, lwd = 2,
       bty = 'n', col = c('royalblue'))


# plot non parametric local regression to for temp increment ------------------------


# span 1 
plot_lowess(temp$year, c(0, diff(temp$temp)), f = 1, ax = temp$year,
            xlab = 'year', ylab = 'temperature', type = 'l',
            main = 'Annual temperature increment')

# add legend
legend('topleft', legend = c('span  100%'), lty = 1, lwd = 2,
       bty = 'n', col = c('violetred'))


# span .5
plot_lowess(temp$year, c(0, diff(temp$temp)), f = .5, ax = temp$year,
            xlab = 'year', ylab = 'temperature', type = 'l',
            main = 'Annual temperature increment',
            add = 0, fit.col = 'royalblue')

# add legend
legend('topleft', legend = c('span  50%'), lty = 1, lwd = 2,
       bty = 'n', col = c('royalblue'))


# join climate data by key = 'year' and compute overall cor -----------------------

X <- merge(rain, temp, by = 'year') # join data in a dataset by year 
cor(X) # compute and print correlation
xtable::xtable(cor(X), type = 'latex') # produce latex format table


# produce pairs plot between all the climate variables  ---------------------------

par(mar = rep(0, 4)) # set margins
pairs(X, col = 'grey60') # shows pairs plot to give an idea of the relationship



# ---------------------------------------------------------------------- #
# ---------------------- climate modelling section --------------------- #
# ---------------------------------------------------------------------- #

# function to numerically integrate ODE 

# model temperature in time with gamma parameter
tempchange <- function(t, state, parms) {
  
  with(as.list(c(state,parms)), { # inner R contenxt
    
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
  out <- ode(y = cinit, times = t, func = tempchange, parms = as.list(parms))
  
  # Filter data that contains time points where data is available
  outdf=data.frame(out) # convert into data frame
  outdf=outdf[outdf$time %in% X$time,] # select time point of interest
  ssqres <- outdf$temp - X$temp # Evaluate predicted vs experimental residual
  
  return(ssqres) # return predicted vs experimental residual
}

# ------------------------------------------------------------------------------------
# fitting LM algorithm to the residual sum of square returned by the function ssq_temp

fitval_temp <- nls.lm(par = parms, fn = ssq_temp) # LM algo

fit_pe <- fitval_temp$par # extract parameters point estimate

# extract parameters standard errors and compute 95% confidence interval ---------------
se <- summary(fitval_temp)$coefficients[1, 'Std. Error']
parest_temp <- c(fitval_temp$par - qnorm(.975) * se, # lower bound
                 fitval_temp$par,                    # point estimate 
                 fitval_temp$par + qnorm(.975) * se) # upper bound



# estimate scenarios with point estimate and confidence interval
out_temp <- ode(y=yini,times=t,func=tempchange,parms=as.list(parest_temp[2])) # point estimate 
out_temp_l <- ode(y=yini,times=t,func=tempchange,parms=as.list(parest_temp[1])) # lower bound
out_temp_u <- ode(y=yini,times=t,func=tempchange,parms=as.list(parest_temp[3])) # upper bound


# plot ODE integration + LM algo result -----------------------------------------------------

ylim <- range(X$temp) # compute y range for plot
plot(out_temp, ylim = ylim, main = 'Annual temperature ODE evolution',
     xlab = 'year', ylab = 'temperature', lwd = 2, col = 'violetred')

# add observed data 
points(X$year - min(X$year), X$temp, col = 'grey70', type = 'l') 

# super impose confidence region
polygon(c(out_temp_l[,1], rev(out_temp_u[,1])),
        c(out_temp_l[,2], rev(out_temp_u[,2])), 
        border = NA, col = ggplot2::alpha('red', .1)
)

# add legend
legend('topleft', col = c('violetred',ggplot2::alpha('red', .1)),
       lty = 1, lwd = c(2, 10), 
       legend = c('ODE estimation', '95% CI'), bty = 'n')




# propose out of sample horizon for 50 years

t_long <- 0:170 #  1901 : 2021 : 2071
X <- merge(X, data.frame(time = t_long), all.y = TRUE) # join data with future years



# numerically integrate ODE with best LM estiamte parameter (95 CI) and longer time
out_temp_lb = ode(y=yini,times=t_long,func=tempchange,parms=as.list(parest_temp[1])) # lower 
out_temp_pe = ode(y=yini,times=t_long,func=tempchange,parms=as.list(parest_temp[2])) # point
out_temp_ub = ode(y=yini,times=t_long,func=tempchange,parms=as.list(parest_temp[3])) # upper

# join all the scenario results by key 'time' in a unique data frame
out_temp_par <- merge(merge(out_temp_lb, out_temp_pe, by = 'time'), out_temp_ub, by = 'time')
names(out_temp_par) <- c('time', 'lb', 'pe', 'ub') # assign name to data set 

# subset only out of sample year, those to be predicted (>2021)
out_temp_par_only <- out_temp_par[out_temp_par$time > max(out_temp[,'time']),]

# plot result
par(mfrow = c(1,1), las = 2) # graphics settings
ylim <- range(out_temp_par[,-1], X$temp, na.rm = TRUE) # ylim range

# plot observed
plot(t_long, X$temp, type = 'l', ylim = ylim, xaxt = 'n', 
     xlab = 'year', ylab = 'temperature', col = 'grey70')

# add point estimate
lines(out_temp_par$time, out_temp_par$pe, col = 'red')

# superimpose confidence region
polygon(c(out_temp_par_only$time, rev(out_temp_par_only$time)),
        c(out_temp_par_only$lb, rev(out_temp_par_only$ub)),
        col = ggplot2::alpha('red', .2), border = NA)

# set axis
axis(1, at = X$time, labels = X$time + min(X$year, na.rm = 1), tick = FALSE)

# add legend
legend('topleft', col = c('violetred',ggplot2::alpha('red', .1)),
       lty = 1, lwd = c(2, 10), 
       legend = c('ODE estimation', '95% CI'), bty = 'n')


# ---------------------------------------------------------------------- #
# -------------------- sensitivity analysis ---------------------------- #
# ---------------------------------------------------------------------- #


# function to numerically integrate ODE + sensitivity equation
tempsens <- function(t, state, parms) {
  
  with(as.list(c(state,parms)), {
    
    gamma <- as.numeric(parms["gamma"])
    RHS1 <- gamma * temp 
    RHS2 <- gamma * s + temp # derived sensitivity equation
    
    return(list(c(RHS2, RHS2)))
  })
}


# loop to explore different scenario: edis 'qs' if you want to explore other scen
qs <- c(5) # number of standard deviation to shift the parameters estimate
f <- function(x) log(x) # choose plot scale: log
 
# loop over proposet
for(q in qs){
  
  sens <- c(yini, s = 0) # set initial condition + 0 for the sensitivity equation
  pe0 <- parest_temp[2]  # set point estimate
  pe1 <- pe0 + q * se    # set upper bound 
  pe2 <- pe0 - q * se    # set lower bound

  # numerically integrate ODE for all the scenario proosed (lb, pe, ub)
  out_temp_sens1 <- ode(y=sens,times=t,func=tempsens,parms=list(gamma = pe0)) # point estimate
  out_temp_sens2 <- ode(y=sens,times=t,func=tempsens,parms=list(gamma = pe1)) # upper bound
  out_temp_sens3 <- ode(y=sens,times=t,func=tempsens,parms=list(gamma = pe2)) # lower bound

  ylim <- f(range(out_temp_sens1[,2], out_temp_sens2[,2], out_temp_sens3[,2])) # y plot range

  # plot time vs temp for all the three chosen scenario (lb, pe, ub)
  plot(out_temp_sens1[,1],  f(out_temp_sens1[,2]), type = 'l', xlab = 'time',
       ylab = '', ylim = ylim)
  lines(out_temp_sens2[,1], f(out_temp_sens2[,2]), type = 'l', col = 'red')
  lines(out_temp_sens3[,1], f(out_temp_sens3[,2]), type = 'l', col = 'blue')
  
  ylim <- f(range(out_temp_sens1[,3], out_temp_sens2[,3], out_temp_sens3[,3]) + 1) # y plot range

  # plot time vs sensitivityf for all the three chosen scenario (lb, pe, ub)
  plot(out_temp_sens1[,1],  f(out_temp_sens1[,3]),
       ylab = '', 
       xlab = 'time', type = 'l', ylim = ylim)
  lines(out_temp_sens2[,1], f(out_temp_sens2[,3]), col  = 'red')
  lines(out_temp_sens3[,1], f(out_temp_sens3[,3]), col  = 'blue')
}




# ---------------------------------------------------------------------- #
# -------------- propose model for rain temperature relationship ------- #
# ---------------------------------------------------------------------- #


rain_temp <- X[complete.cases(X),] # subset X where there is no NA value

# propose a linear model for rain temp relationship plus non linear spline term
# spline K = 2 has empirially shown to be the best for performance and parsimony
# temp_sd has shown to be not significant, hence we will not include it in the model

fit_rain <- lm(rain ~ splines::ns(temp, 2), data = rain_temp) # fit model
BIC(fit_rain) # compute BIC
summary(fit_rain) # print model summary

# inspect fitted vs residuals diagnostic:
plot(fitted(fit_rain), residuals(fit_rain), main = 'OLS',
     xlab = 'fitted', ylab = 'residuals', col = 'grey70')
legend('topleft', legend = paste('BIC', round(BIC(fit_rain),2)), bty = 'n')


# to remedy, iterative WLS for 5 times

fit_rain_wi <- fit_rain # copy model object
for(i in 1:5){ # iterate 5 times
  wi <- 1 / residuals(fit_rain_wi)**2 # extract residuals from last model and compute weights
  fit_rain_wi <- lm(rain ~ splines::ns(temp, 2), data = rain_temp, weights = wi) # fit model with weights
}

# plot WLS model results
plot(fitted(fit_rain_wi), residuals(fit_rain_wi), main = 'WLS',
     xlab = 'fitted', ylab = 'residuals', col = 'grey70')
legend('topleft', legend = paste('BIC', round(BIC(fit_rain_wi),2)), bty = 'n')

# combina and compare OLS and WLS coefficients standard error
ols_wls_se <-
  cbind(
    summary(fit_rain)$coefficients[,2],
    summary(fit_rain_wi)$coefficients[,2]
  )
# xtable::xtable(ols_wls_se, type = 'latex', digits = 8) # produce table for latex format





# ---------------------------------------------------------------------- #
# -------------- sensitivity analysis for rain temp model -------------- #
# ---------------------------------------------------------------------- #


# starting from temperature lower and upper bound we compute either lower and upper bound
# by doing so we include in the model the DOUBLE UNCERTAINTY source and store the resulting
# lower, upper bound and point estimate 

# predict rain from temp
pred_type <- 'prediction' # prediction type will produce confidence interval
rain_pred_pe <- predict(fit_rain_wi, newdata = data.frame(temp = out_temp_par$pe), interval = pred_type)[,'fit']
rain_pred_lb <- predict(fit_rain_wi, newdata = data.frame(temp = out_temp_par$lb), interval = pred_type)[,'lwr']
rain_pred_ub <- predict(fit_rain_wi, newdata = data.frame(temp = out_temp_par$ub), interval = pred_type)[,'upr']

# sensitivity
# try to impute different values
q <- 5 # number of standard error to shift the parameter estimate

fit_rain_wi_edit <- fit_rain_wi # copy model object 
fit_rain_wi_edit$coefficients[2] <-  -8.962673 + q *  1.0267 # shift parameters -> upper
fit_rain_wi_edit$coefficients[3] <-  -8.6207 + q *  0.3338 # shift parameters -> upper

# repeat prediction with lowewr shifted model
pred_type <- 'prediction'
rain_pred_pe_e <- predict(fit_rain_wi_edit, newdata = data.frame(temp = out_temp_par$pe), interval = pred_type)[,'fit']
rain_pred_lb_e <- predict(fit_rain_wi_edit, newdata = data.frame(temp = out_temp_par$lb), interval = pred_type)[,'lwr']
rain_pred_ub_e <- predict(fit_rain_wi_edit, newdata = data.frame(temp = out_temp_par$ub), interval = pred_type)[,'upr']

# edit model
fit_rain_wi_edit$coefficients[2] <-  -8.962673 - q *  1.0267 # shift parameters -> lower
fit_rain_wi_edit$coefficients[3] <-  -8.6207 - q *  0.3338 # shift parameters -> lower

# repeat prediction with upper shifted model
pred_type <- 'prediction'
rain_pred_pe_f <- predict(fit_rain_wi_edit, newdata = data.frame(temp = out_temp_par$pe), interval = pred_type)[,'fit']
rain_pred_lb_f <- predict(fit_rain_wi_edit, newdata = data.frame(temp = out_temp_par$lb), interval = pred_type)[,'lwr']
rain_pred_ub_f <- predict(fit_rain_wi_edit, newdata = data.frame(temp = out_temp_par$ub), interval = pred_type)[,'upr']


# plot the result
par(las = 2) # graph param
plot(rain_temp$time, rain_temp$rain, type = 'l', xlim = range(out_temp_par$time),
     xaxt = 'n', col = 'grey70', xlab = 'time', ylab = 'rain')

# add axis
axis(1, at = out_temp_par$time, labels = out_temp_par$time + min(X$year, na.rm = TRUE), tick = FALSE)
title('Rain vs Time') # add title

# add point estimate
lines(out_temp_par$time, rain_pred_pe, lwd = 2, col = 'red')

# subset out of sample points
idx <- out_temp_par$time > max(rain_temp$time)

# add confidence region for out of sample points
polygon(c(out_temp_par$time[idx], rev(out_temp_par$time[idx])),
        c(rain_pred_lb[idx], rev(rain_pred_ub[idx])),
        border = NA, col = ggplot2::alpha('red', .2))

# add confidence region for other scenario + se
lines(out_temp_par$time, rain_pred_pe_e, lwd = 2, col = 'blue')
polygon(c(out_temp_par$time[idx], rev(out_temp_par$time[idx])),
        c(rain_pred_lb_e[idx], rev(rain_pred_ub_e[idx])),
        border = NA, col = ggplot2::alpha('blue', .1))

# add confidence region for other scenario - se
lines(out_temp_par$time, rain_pred_pe_f, lwd = 2, col = 'green')
polygon(c(out_temp_par$time[idx], rev(out_temp_par$time[idx])),
        c(rain_pred_lb_f[idx], rev(rain_pred_ub_f[idx])),
        border = NA, col = ggplot2::alpha('green', .1))

# add legend
legend('topright', col = c('red','blue','green'), bty = 'n', lty = 1, lwd = 2,
       legend = c('Estimated pars',
                  paste('est +', round(q, 2), 'SE'),
                  paste('est -', round(q, 2), 'SE')
       )
)





# --------------------------------------------------------------------------- #
# ----------------------- Species count Modelling --------------------------- #                        #
# --------------------------------------------------------------------------- #



# count data have been splitted to be stored on github
# load and row bind the different sets in a unique data frame
dat0 <- rbind(
  read.csv("../data/occurrence1.txt", row.names=1),
  read.csv("../data/occurrence2.txt", row.names=1),
  read.csv("../data/occurrence3.txt", row.names=1),
  read.csv("../data/occurrence4.txt", row.names=1)
)


dat <- dat0 # backup data 

# remove empty species
dat <- dat[ - which(dat$family == ''), ] 

# sanity check on year: keep only 1900+
dat <- dat[dat$year > 1900,]
years <- cbind( year = min(dat$year) : max(dat$year) ) # create all range year sequence

# create temp dat to compute annual mean data frame
temp <- dat # copy object

# produce a COMPOUND KEY which will be use within the tapply function
# we want to distinguish each species and each genes per year
temp$key <- paste0( dat$year, '_', dat$family, '.', dat$genus) # create KEY 

# compute annual mean for year : species : genes
counts <- tapply(temp$individualCount, temp$key, mean) 

# store results in a data frame
Y <- data.frame(
  year = as.numeric(substr(names(counts), 1, 4)), # extract year from key
  species = unlist(lapply(strsplit(names(counts), '_'), '[[', 2)), # extract species.genes from key
  individualCount = counts # attach yearly count
)

# count

# subset data by removing NA cases and start from 1975
Y <- Y[complete.cases(Y) & Y$year > 1900,]

fams_tab <- table(Y$species) # extract unique species.genes

# produce barpot of the number of families
par(las = 2, mar = c(12, 5, 2, 5)) # graph param
barplot(fams_tab[order(fams_tab, decreasing = TRUE)], border = NA,
        col = ggplot2::alpha(viridis(length(fams_tab), direction = -1), .4) )
abline(h = 30, col = 'royalblue', lwd = 1, lty = 2) # set horizontal line at 30

# subset observation with n > 30
obs_30 <- names(fams_tab[which(fams_tab > 30)])
Y <- Y[Y$species %in% obs_30,]


# pivot table: convert long table into wide table
# where the first column: year 
# and one column per each species.genes combination

Y_wide <- data.frame(
  tidyr::pivot_wider(
    Y,
    names_from = 'species',
    values_from = 'individualCount'
  )
)

# quick visualization
matplot(Y_wide[,1], Y_wide[,-1], type = 'l')

# remove extreme observations -----------------------------------------------

# compute stats: mean,  sd
Y_mean <- mean(Y$individualCount, na.rm = TRUE)
Y_sd <- sd(Y$individualCount, na.rm = TRUE)

# identify those above/below the threshold
q <- qnorm(.99) # quantile to look at 
out_idx <-
  which(
    Y$individualCount < (Y_mean  - q * Y_sd) | 
      Y$individualCount > (Y_mean  + q * Y_sd)
  )

# identify species to be removed
species_drop <-
  sort(
    unique(
      Y$species[Y$species %in% unique(Y$species[out_idx])]
    )
  )

# visualize: first plot to be remove, second plot those left
Y_wide_col <- which(names(Y_wide) %in% species_drop)
matplot(Y_wide[,1], Y_wide[,   Y_wide_col], type = 'l')
matplot(Y_wide[,1], Y_wide[, - Y_wide_col], type = 'l', ylim = c(0, 1300), lwd = .4)


# remove species from both Y_wide and Y
Y_wide <- Y_wide[ - Y_wide_col]
Y <- Y[ - which(Y$species %in% species_drop),]



# provide a quick visualization ----------------------------------------------

# visulize with mean and median
wide_mean <- apply(Y_wide[,-1], 1, function(x) mean(x, na.rm = TRUE))
wide_median <- apply(Y_wide[,-1], 1, function(x) median(x, na.rm = TRUE))

# compute local weighted regression non parametric fit for both mean and median
wide_mean_fit <- fitted(loess(wide_mean ~ Y_wide[,1]))
wide_median_fit <- fitted(loess(wide_median ~ Y_wide[,1]))

# viz with median and mean non parametric fit 
par(las = 1, mar = c(2, 5, 0, 5))
layout(matrix(c(1,1,1,1,2,2), ncol = 2, byrow = 1))
matplot(Y_wide[,1], Y_wide[,-1], type = 'l', lwd = .2, lty = 1,
        xlab = 'year', ylab = 'individual count', xaxt = 'n')
lines(Y_wide[,1], wide_mean_fit, lwd = 2, col = 'red')
lines(Y_wide[,1], wide_median_fit, lwd = 2, col = 'blue')
matplot(Y_wide[,1], Y_wide[,-1], type = 'l', lwd = .2, ylim = c(0, 100),
        xlab = 'year', ylab = 'individual count', xaxt = 'n', lty = 1)
lines(Y_wide[,1], wide_mean_fit, lwd = 2, col = 'red')
lines(Y_wide[,1], wide_median_fit, lwd = 2, col = 'blue')



# load climate data into the environment
X <- read.csv('../data/climate_pred.csv')[,-1]

# combine climate data and species count by year ---------------------------------------------

Y <- merge(Y, X, by = 'year', all.x = TRUE) # join data by year
Y <- Y[order(Y$species), ] # sort data by species
Y$individualCount <- round(
  as.numeric(Y$individualCount)
) # cast into integere


# inspect individual count distribution: we can try a Pois distr
barplot(table(Y$individualCount) / nrow(Y), border = NA, col = viridis(600),
        xlab = 'individual count', ylab = 'frequencies')


# data preprocessing: shift data so that it is centered on zero ---------------------------
cyear <- 1998
Y$syear <- Y$year - cyear

# standardize rain: mean 0 variance 1 and store result to be applied to out out sample obs
m_rain <- mean(Y$rain_obs, na.rm = TRUE) # compute mean
sd_rain <- sd(Y$rain_obs, na.rm = TRUE) # compute sd
Y$srain <- (Y$rain_obs - m_rain)/sd_rain # (x-m)/s


# standardize temp: mean 0 variance 1 and store result to be applied to out out sample obs
m_temp <- mean(Y$temp_obs, na.rm = TRUE) # compute mean
sd_temp <- sd(Y$temp_obs, na.rm = TRUE) # compute sd 
Y$stemp <- (Y$temp_obs - m_temp)/sd_temp # (x-m)/s


# empirical trial to choose the optimum spline degree --------------------------
# here we fit the model where the reponse is the individual count for each 
# species-genes pairs and the predictors are the climate data (rain, temp)
# and some interactions, those which have shown to be significant

# # define degree
# bics <- rep(NA, 7)
# for(i in 1:7){
#   fit_year_temp <- lme4::g\textilmer(
#     # this model returns the best AIC among those explored
#     individualCount ~ ns(syear, i) * stemp +
#       ns(syear, i) * srain + (1 + stemp| species), 
#     family = 'poisson', data = Y
#   )
#   bics[i] <- BIC(fit_year_temp)
# }
# plot(1:7, bics, type = 'b', pch = 20)

q_temp <- 5 # choose spline degree 

# fit model
fit_year_temp <- lme4::glmer(
  # this model returns the best AIC among those explored
  individualCount ~ ns(syear, q_temp) * stemp +
    ns(syear, q_temp) * srain + (1 + stemp| species), 
  family = 'poisson', data = Y
)

# store model summary
model_res <- summary(fit_year_temp)
xtable::xtable(model_res$coefficients, type = 'latex') # latex table format


# model diagnostic: as it is possion we expect most of the fitted to be low values
# and the residuals to be randomly distributed across zero. we also expect
# the variance to increas with the mean, as per poisson assumption

plot(fitted(fit_year_temp), residuals(fit_year_temp),
     xlab = 'fitted values', ylab = 'residuals', col = 'grey40', lwd = .25)

# extract random effect from the model
fit_ranef <- lme4::ranef(fit_year_temp, condVar = TRUE)
fit_ranef_df <- fit_ranef$species

# check qqplot for each of the random effect
P <- ncol(fit_ranef_df)
mux <- sdx <- rep(NA, P)
for(i in 1:P){
  
  # normality check
  qqnorm(fit_ranef_df[,i], main = names(fit_ranef_df)[i])
  qqline(fit_ranef_df[,i], col = 'red', lwd = 1.5)
  
}

# extract model terms contribution to variance
sigma_eps <- 1.000
sigma_0 <- 1.3586  
sigma_temp <- 0.2597 

# combine the variances into a vector
vars <- c(sigma_0, sigma_temp, sigma_eps)
names(vars) <- c('intercept','temperature','residuals') # assign names


# plot percentage contirbution of each term to total variance
barplot(vars / sum(vars), border = NA, col = 'grey80',
        main = 'Random effect variance contribution')


# plot confidence interval for each of the random effect term (intercept + slope)
dotplot_ranef <- lattice::dotplot(fit_ranef)
dotplot_ranef$species



# extract fitted values from the model and attach to the original data set
Y$fitted[complete.cases(Y)] <- fitted(fit_year_temp)
spec_tab <- table(Y$species) # compute species numerousness
rmse <- rep(NA, length(spec_tab)) # initialize vector to store rmse (percentage)


##################### RMSE hasn't been included in final report ############################

# for each species compute median rmse percentage and store into vector
for(i in 1:length(spec_tab)){
  
  # subset  
  ispec <- names(spec_tab[i]) # select species
  iy <- Y[Y$species == names(spec_tab[i]),] # subset data
  
  # compute rmse
  rmse[i] <- median(
    ( (iy$individualCount - iy$fitted) / iy$individualCount )^2,
    na.rm = TRUE
  )
}
rmse <- sqrt(rmse) # compute square root
names(rmse) <- names(spec_tab) # assign names to vector
order_idx <- order(rmse, decreasing = 1) # compute decreasing order indexes

# plot the results
par(las = 2, mar = c(14, 4, 2, 2), mfrow = c(1,1))
barplot(rmse[order_idx], main = 'median RMSE', border = NA)

#############################################################################################

# apply pre processing to X so that it matches the model data frame
X$year <- min(X$year, na.rm = TRUE) + X$time # compute long years
X$syear <- X$year - cyear # center year
X$srain_lb <- (X$rain_pred_lb - m_rain) / sd_rain # std rain lb
X$srain_ub <- (X$rain_pred_ub - m_rain) / sd_rain # std rain ub
X$stemp_lb <- (X$temp_pred_lb - m_temp) / sd_temp # std temp lb
X$stemp_ub <- (X$temp_pred_ub - m_temp) / sd_temp # std temp ub
X$srain <- (X$srain_ub - X$srain_lb)/2 + X$srain_lb # compute rain point estimate
X$stemp <- (X$stemp_ub - X$stemp_lb)/2 + X$stemp_lb # compute temp point estimate
X <- X[X$year >= 1975,] # subset data to those higher than 1975


# extract coefficients SE from model estimation
se <-
  sqrt(
    diag(
      as.matrix(
        vcov(
          fit_year_temp
        )
      )
    )
  )

# compute model 95 confidence interval for model coefficients
q <- qnorm(.975) # compute normal quantile
tab <- data.frame( # store in a data frame
  lb = lme4::fixef(fit_year_temp) - q * se,
  pe = lme4::fixef(fit_year_temp),
  ub = lme4::fixef(fit_year_temp) + q * se
)

# create both lb and ub copy model
fit_year_temp_lb <- fit_year_temp_ub <- fit_year_temp
fit_year_temp_lb@beta <- tab$lb # assign lower bounded coefficients
fit_year_temp_ub@beta <- tab$ub # assign upper bounded coefficients

# extract species names
spec <- names(rmse)


par(mfrow = c(3,3), las = 2, mar = c(3,3,2,0))

# for each of the species - genes involved
for(i in 1:length(spec)){
  
  # predict upper and lower bound according to species
  pred_lb <-  predict(
    fit_year_temp_lb,
    data.frame(
      syear = X$syear,
      stemp = X$stemp_lb,
      srain = X$srain_lb,
      species = spec[i]
    ),
    type = 'response'
  )
  pred_ub <-  predict(
    fit_year_temp_ub,
    data.frame(
      syear = X$syear,
      stemp = X$stemp_ub,
      srain = X$srain_ub,
      species = spec[i]
    ),
    type = 'response'
  )
  
  # compute point estimate
  pred_pe <- (pred_ub - pred_lb)/2 + pred_lb
  
  
  # ylim range values
  ylim <- c(0, 
            max(
              Y$individualCount[Y$species == spec[i]],
              pred_lb, pred_ub
            )
  )
  
  # plot observed data
  plot(
    X$year,
    c(
      Y$individualCount[Y$species == spec[i]],
      rep(NA, nrow(X) - sum(Y$species == spec[i]))
    ),
    type = 'l',
    ylim = ylim,
    xlab = '',
    ylab = paste(spec[i],'individual count'),
    main = paste0(i,': ',spec[i]),
    xaxt = 'n',
  )
  
  
  # plot axis value at bottom line
  if(i %% 18 > 15 | i %% 18 == 0) axis(1, at = X$year, labels = X$year, tick = FALSE)
  
  # superimpose estimated confidence interval
  polygon(
    c(X$year, rev(X$year)),
    c(pred_lb, rev(pred_ub)),
    border = NA,
    col = ggplot2::alpha(2, .2)
  )
  
  
  # superimpose point estimate
  lines(
    X$year,
    pred_pe,
    lwd = 1.5,
    col = 'red')
  
  # indentify 1 cases
  zero_case <- which(round(pred_lb) == 1)
  zero_case <- zero_case[zero_case > 2021-1975]

  # add vertical line to indicate 1 individual hit
  if(length(zero_case) > 0) {
    abline(v = X$year[zero_case[1]], lty = 2)
    legend('topright', col = 'grey', lty = 2, legend = X$year[zero_case[1]], bty = 'n')
    hit_one[i] <- X$year[zero_case[1]] # store result
  }
}


# compute stas
mean(hit_one)
median(hit_one)