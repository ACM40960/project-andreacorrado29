
# --------------------------------------------------------------------------- #
#                               math project final                            #
# --------------------------------------------------------------------------- #

# diff_eq_cov.R
# .rs.restartR()
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

rain <- read.csv('pr_1901_2020_ZAF.csv') # load data
names(rain)[1] <- 'rain' # assign name
rain_mean <- tapply(rain$rain, rain$Year, mean)
rain_sd <- tapply(rain$rain, rain$Year, sd)
rain <- data.frame(year = as.numeric(names(rain_mean)), rain = rain_mean, rain_sd = rain_sd)
plot_lowess(rain$year, rain$rain, f = .75, main = 'yearly rain', ax = rain$year)


# --------------------------------------------------------------------------- #
#                                temperature                                  #
# --------------------------------------------------------------------------- #

temp <- read.csv('tas_1901_2020_ZAF.csv')
names(temp)[1] <- 'temp'
temp_mean <-tapply(temp$temp, temp$Year, mean)
temp_sd <-tapply(temp$temp, temp$Year, sd)
temp <- data.frame(year = names(temp_mean), temp = temp_mean, temp_sd = temp_sd)
plot_lowess(temp$year, temp$temp, f = .75, main = 'yearly temp', ax = temp$year)

# --------------------------------------------------------------------------- #
#                               co2 emissions                                 #
# --------------------------------------------------------------------------- #


co2 <- read.csv("co2_raw.txt")
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

if(!('dat0' %in% ls())){
  dat0 <- read.delim("~/Dropbox/ucd_MATH_PROJ/DROUGHT/occurrence.txt") # load data
}

dat <- dat0
dat <- dat[ - which(dat$family == ''), ] # remove empty species

# sanity check on year
sort(unique(dat$year))
dat <- dat[dat$year > 1900,]
years <- cbind( year = min(dat$year) : max(dat$year) )

# pre process data
Y <- data.frame()
min.row <- 45
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

# need to remove Scolopacidae as the values are too high compared to 
# the other species and we are not going to standardize or take relative change
# as we are interested in the absolute number
# Y <- Y[ - which(Y$species == 'Scolopacidae'), ]

# wrong values
Y <- Y[- which(Y$species == 'Scolopacidae' & Y$year == 1975), ] 

# pivot table
Y_wide <- data.frame(
  tidyr::pivot_wider(Y,
                     names_from = species,
                     values_from = individualCount
                     )
  )

# plot families
par(las = 2)
matplot(Y_wide[,1], Y_wide[,-1], type = 'l', xlab = 'year', ylab = 'individualCount', xaxt = 'n')
axis(1, at = Y_wide[,1], labels = Y_wide[,1])
title('Species count vs time')

# --------------------------------------------------------------------------- #
#                            model of interest                                #
# --------------------------------------------------------------------------- #

# combine data ----------------------------------------------------------------
Y <- merge(Y, X, by = 'year', all.x = TRUE)
Y <- Y[order(Y$species), ]

# center data -----------------------------------------------------------------
scale <- TRUE
Y$syear <- scale(Y$year, scale = scale)
Y$stemp <- scale(Y$temp, scale = scale)
Y$srain <- scale(Y$rain, scale = scale)
Y$individualCount <- round(as.numeric(Y$individualCount))
Y <- Y[!is.na(Y$temp),]

# fit model -------------------------------------------------------------------
fit <-
  lme4::glmer(
    individualCount ~ stemp * year + srain * year + (1 | species), 
    family = 'poisson', data = Y
)
summary(fit) 


# visualize result ------------------------------------------------------------

Y$fitted <- fitted(fit)
Y_wide <- data.frame(
  tidyr::pivot_wider(Y,
                     names_from = species,
                     values_from = c(individualCount, fitted)
  )
)

 
matplot(Y_wide$year, Y_wide[,grep('individualCount', names(Y_wide))], type = 'p',
        xlab = 'years', ylab = 'individualCount', xaxt = 'n')
matlines(Y_wide$year, Y_wide[,grep('fitted', names(Y_wide))], type = 'l')
axis(1, at = Y_wide[,1], labels = Y_wide[,1])



#
# we now want to focus on one of the species, we take the one for which
# the valeus are not zero yet

Y1 <- Y[Y$species == 'Charadriidae', - which(names(Y) == 'species')]

formula(fit) # individualCount ~ stemp * srain * syear - stemp * srain 

fitg <- 
  glm(individualCount ~ stemp * year + srain * year, 
      data = Y1, family = 'poisson')
summary(fitg)

# -----------------------------------------------------------------------------

Y1$fitted <- fitted(fitg)


plot(Y1$year, Y1$individualCount, ylim = c(0, max(Y1$individualCount)), xaxt = 'n',
     xlab = 'year', ylab = 'individualCount')
title('Charadriidae in time')
lines(Y1$year, Y1$fitted, type = 'l', col = 'red')
axis(1, at = Y_wide[,1], labels = Y_wide[,1])

years_plus <- 1975:2030
t_long <- 0:length(years_plus) # times long 
Y1$time <- 0 : ( nrow(Y1) - 1)

# --------------------------------------------------------------------------- #
#                      model rain data with differential eq                   #
# --------------------------------------------------------------------------- #

# rain
rainchange <- function(t,c,parms){# rate constant passed through a list called parms
  
  alpha <- parms$alpha
  
  # derivatives dr/dt are computed below
  r <- rep(0,length(c))
  r[1] <- alpha * c["rain"]
  
  return(list(r))
}

warning('write as linear system ?')
rainchange <- function(t, state, parms) {

    with(as.list(c(state,parms)), {
      
      alpha <- as.numeric(parms["alpha"]) 
      beta <- as.numeric(parms["beta"]) 
      gamma <- as.numeric(parms["gamma"]) 
    
      dR = alpha * rain + beta * temp
      #dR = (alpha + beta * temp) * rain
      dT = gamma * temp
    
      return(list(c(dR,dT)))
  })
}

t <- Y1$time
yini <- c(rain = Y1$rain[1], temp = Y1$temp[1])
parms <- c(alpha = -1e-2, beta = 1e-2, gamma = 1e-3)
out <- deSolve::ode(y = yini, time = t, func = rainchange, parms = as.list(parms))
out

# function that calculates residual sum of squares
ssq_rain <- function(parms){
  
  #browser()
  
  # inital concentration
  cinit <- c(rain = Y1$rain[1], temp = Y1$temp[1])
  # time points for which conc is reported
  t <- sort( unique( c(seq(0,5,0.1),Y1$time) ))
  
  # solve ODE for a given set of parameters
  out=ode(y=cinit,times=t,func=rainchange,parms=as.list(parms))
  
  # Filter data that contains time points where data is available
  outdf=data.frame(out)
  outdf=outdf[outdf$time %in% Y1$time,]
  # Evaluate predicted vs experimental residual
  ssqres <- c(outdf$rain-Y1$rain, outdf$temp - Y1$temp)
  
  return(ssqres) # return predicted vs experimental residual
}

# parameter fitting using levenberg marquart algorithm
parms_rain <- c(alpha = -1e-2, beta = 1e-2, gamma = 1e-3)

# fitting
fitval_rain <- nls.lm(par = parms_rain, fn = ssq_rain)
summary(fitval_rain)

parest_rain <- fitval_rain$par # parest
cinit <- c(rain = Y1$rain[1], temp = Y1$temp[1])
out_rain = ode(y=cinit,times=t_long,func=rainchange,parms=as.list(parest_rain))
out_rain <- data.frame(out_rain)

par(mfrow = c(1,2))
plot(Y1$time, Y1$rain, type = 'b')
points(out_rain$time, out_rain$rain, col = 'red')
plot(Y1$time, Y1$temp, type = 'b')
points(out_rain$time, out_rain$temp, col = 'red')


# --------------------------------------------------------------------------- #
#                      model temp data with differential eq                   #
# --------------------------------------------------------------------------- #


# rain
tempchange <- function(t,c,parms){# rate constant passed through a list called parms
  
  beta <- parms$beta
  
  # derivatives dr/dt are computed below
  r <- rep(0,length(c))
  r[1] <- beta * c["temp"]
  
  return(list(r))
}

# function that calculates residual sum of squares
ssq_temp <- function(parms){
  #browser()
  
  # inital concentration
  cinit <- c(temp = Y1$temp[1])
  t=c(seq(0,5,1),Y1$time)
  t=sort(unique(t))
  
  # solve ODE for a given set of parameters
  out=ode(y=cinit,times=t,func=tempchange,parms=as.list(parms))
  
  # Filter data that contains time points where data is available
  outdf=data.frame(out)
  outdf=outdf[outdf$time %in% Y1$time,]
  ssqres=outdf$temp-Y1$temp
  
  return(ssqres) # return predicted vs experimental residual
}

# parameter fitting using levenberg marquart algorithm
cinit <- c(temp = Y1$temp[1])
t <- Y1$time
parms_temp <- c(beta = 0.1)

# fitting
fitval_temp <- nls.lm(par = parms_temp, fn = ssq_temp)
summary(fitval_temp)

parest_temp <- fitval_temp$par # parest
out_temp = ode(y=cinit,times=t_long,func=tempchange,parms=as.list(parest_temp))

# visualize results ------------------------------------------------------------

par(las = 2, mfrow = c(1,2))

plot(out_rain[,1], out_rain[,2], xaxt = 'n', lwd = 2, col = 'royalblue', type = 'l',
     xlab = 'time', ylab = 'rain')
title('Evolution of rain in time')
points(Y1$time, Y1$rain)
axis(1, at = t_long, t_long + 1975)


plot(out_temp[,1], out_temp[,2], xaxt = 'n', lwd = 2, col = 'darkorange2', type = 'l',
     xlab = 'time', ylab = 'temperature')
title('Evolution of temperature in time')
points(Y1$time, Y1$temp)
axis(1, at = t_long, t_long + 1975)



# i want to estimate the uncertainty of these estimation
 
warning('do you really want to specify your matrix in this way?')
S <- diag(c(vcov(fitval_rain), vcov(fitval_temp)))
dof <- nrow(Y1) * 3 #dof

# draw the confidence region
# get points for a circle with radius r
r=sqrt(qf(0.95,2,dof)*2)
theta=seq(0,2*pi,length.out=100)
z=cbind(r*cos(theta),r*sin(theta))
# transform points of circle into points of ellipse using
# svd of covariance matrix
S_svd=svd(S)      # SVD of covariance matrix
xt=t(S_svd$v)%*%diag(sqrt(S_svd$d))%*%t(z) # transform from circle to ellispse
x=t(xt)
# translate the ellipse so that center is the estimated parameter value
parest <- c(alpha = parest_rain, beta = parest_temp)
x=x+matrix(rep(as.numeric(parest),100),nrow=100,byrow=T)

par(mfrow = c(1,1))
plot(x[,1],x[,2],type="l",xlab="alpha",ylab="beta",lwd=2)
points(parest[1], parest[2], pch=20,col="blue",cex=2)



















