 # count_model.R

rm(list = ls()) # clean env
source('functions.R') # load functions  

# laod packages: need to be installed if not available yet
library(ggplot2) 
library(gam)
library(reshape2) 
library(deSolve) 
library(minpack.lm)


# --------------------------------------------------------------------------- #
#                            load count data                                  #
# --------------------------------------------------------------------------- #

dat0 <- rbind(
  read.csv("../data/occurrence1.txt", row.names=1),
  read.csv("../data/occurrence2.txt", row.names=1),
  read.csv("../data/occurrence3.txt", row.names=1),
  read.csv("../data/occurrence4.txt", row.names=1)
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

# remove the NA cases
# compute annual mean
# fit a non parametric local regression to understand the trend
# quick viz
Y_complete <- Y[complete.cases(Y),]
Y_complete$individualCount <- as.numeric(Y_complete$individualCount)

means <- tapply(Y_complete$individualCount, Y_complete$year, mean)
years <- as.numeric(names(means))

# quick visualization
Y_wide <- data.frame(Y_wide)
matplot(Y_wide[,1], Y_wide[,-1], type = 'l', xaxt = 'n', lwd = .5)
lines(years, fitted(loess(means ~ years, span = 1)), lwd = 2, col = 'red')
axis(1, at = Y_wide[,1], labels = Y_wide[,1])
title('Species count vs time')

# --------------------------------------------------------------------------- #
#                            model of interest                                #
# --------------------------------------------------------------------------- #

# stop('before reading climate data')

# load climate data
X <- read.csv('../data/climate_pred.csv')[,-1]

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
par(mfrow = c(1,1))

# center data -----------------------------------------------------------------
cyear <- 1998
Y$syear <- Y$year - cyear

# fit model -------------------------------------------------------------------


# define degree
Q <- 6
bics <- rep(NA, Q)
for(q in 1:Q){
  q_fit_year <- lme4::glmer(
    individualCount ~ ns(syear, q) + (1 | species), 
    family = 'poisson', data = Y
  )
  bics[q] <- BIC(q_fit_year)
}
# balance between model degree and interpretability
q <- 4
cols <- rep('royalblue', Q); cols[q] <- 'red'
plot(1:Q, bics, type = 'b', col = cols, pch = 20)

# fit model
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
# therefore we will try include further predictors
plot(fitted(fit_year), residuals(fit_year))

# further pred: preprocessing
m_rain <- mean(Y$rain_obs, na.rm = TRUE)
sd_rain <- sd(Y$rain_obs, na.rm = TRUE)
Y$srain <- (Y$rain_obs - m_rain)/sd_rain

m_temp <- mean(Y$temp_obs, na.rm = TRUE)
sd_temp <- sd(Y$temp_obs, na.rm = TRUE)
Y$stemp <- (Y$temp_obs - m_temp)/sd_temp


# Y$srain <- scale(Y$rain_obs, scale = TRUE)
# Y$stemp <- scale(Y$temp_obs, scale = TRUE)

Q_temp <- 6
bics <- rep(NA, Q_temp)
for(q_temp in 1:Q_temp){
  q_fit_year_temp <- lme4::glmer(
    # this model returns the best AIC among those explored
    individualCount ~ ns(syear, q_temp) * stemp + ns(syear, q_temp) * srain + (1 | species), 
    family = 'poisson', data = Y
  )
  bics[q_temp] <- BIC(q_fit_year_temp)
}

# balance between model degree and interpretability: again 4 may be a good choice
q_temp <- 4
cols <- rep('royalblue', Q_temp); cols[q_temp] <- 'red'
plot(1:Q_temp, bics, type = 'b', col = cols, pch = 20)

fit_year_temp <- lme4::glmer(
  # this model returns the best AIC among those explored
  individualCount ~ ns(syear, q_temp) * stemp + ns(syear, q_temp) * srain + (1 | species), 
  family = 'poisson', data = Y
)
summary(fit_year_temp)

# less pronounced but still some structure 
par(mfrow = c(1,2))
plot(fitted(fit_year), residuals(fit_year))
plot(fitted(fit_year_temp), residuals(fit_year_temp))

# decreased by 5%
cbind(
  BIC(fit_year),
  BIC(fit_year_temp), # lower BIC
  BIC(fit_year_temp) / BIC(fit_year)
)

# stop('2146 20210722')

# --------------------------------------------------------------------------- #
#                     focus on subset for further analysis                   #
# --------------------------------------------------------------------------- #


# species count subset: take the one with most observations
to_sub <- names(
  sort(table(Y[!is.na(Y$individualCount),'species']), decreasing = TRUE)[1]
)
Y1 <- Y[Y$species == to_sub & complete.cases(Y),] # remove NA obs
dim(Y1)

par(mfrow = c(1,3))
plot(Y1$year, Y1$individualCount, type = 'l')
plot(Y1$year, Y1$srain, type = 'l')
plot(Y1$year, Y1$stemp, type = 'l')

# fit glm
fit_glm <- glm(individualCount ~ ns(syear, q_temp) * stemp + ns(syear, q_temp) * srain,
               data = Y1, family = 'poisson')
summary(fit_glm)
plot(fitted(fit_glm), residuals(fit_glm))


# preprocess X
head(X)
tail(X)

X$year <- X$time + 1901
X$syear <- X$year - cyear
X$srain_lb <- (X$rain_pred_lb - m_rain) / sd_rain
X$srain_ub <- (X$rain_pred_ub - m_rain) / sd_rain
X$stemp_lb <- (X$temp_pred_lb - m_temp) / sd_temp
X$stemp_ub <- (X$temp_pred_ub - m_temp) / sd_temp
X <- X[X$year >= 1975,]


# predict upper and lower bound
pred_lb <- predict(fit_glm, se.fit = TRUE,
                   newdata = data.frame(syear = X$syear,
                                        stemp = X$stemp_lb,
                                        srain = X$srain_lb)
)
pred_ub <- predict(fit_glm, se.fit = TRUE,
                   newdata = data.frame(syear = X$syear,
                                        stemp = X$stemp_ub,
                                        srain = X$srain_ub)
)
# compute lower and upper bound
pred_lb <- pred_lb$fit - qnorm(.975) * pred_lb$se.fit
pred_ub <- pred_ub$fit + qnorm(.975) * pred_ub$se.fit

# apply g transformation
g_pred_lb <- exp(pred_lb)
g_pred_ub <- exp(pred_ub)


newX <- merge(cbind(year = X$year, g_pred_lb, g_pred_ub), 
              Y1[,c('year','individualCount')], 
              by = 'year',
              all.x = TRUE)

par(mfrow = c(1,1))
plot(newX$year, newX$individualCount, type = 'l')
polygon(c(newX$year, rev(newX$year)), c(newX$g_pred_lb, rev(newX$g_pred_ub)),
        border = NA, col = ggplot2::alpha('red', .2))

stop('arrivede here 20210729 2140')

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



















