# practical example

# x pre processing
# X$year <- min(X$year, na.rm = TRUE) + X$time
# X$syear <- X$year - cyear
# X$srain_lb <- (X$rain_pred_lb - m_rain) / sd_rain
# X$srain_ub <- (X$rain_pred_ub - m_rain) / sd_rain
# X$stemp_lb <- (X$temp_pred_lb - m_temp) / sd_temp
# X$stemp_ub <- (X$temp_pred_ub - m_temp) / sd_temp
# X$srain <- (X$srain_ub - X$srain_lb)/2 + X$srain_lb
# X$stemp <- (X$stemp_ub - X$stemp_lb)/2 + X$stemp_lb
# X <- X[X$year >= 1975,]
# head(X)
# tail(X)


# extract se
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

# compute model 95 confidence interval
q <- qnorm(.975)
tab <- data.frame(
  lb = fixef(fit_year_temp) - q * se,
  pe = fixef(fit_year_temp),
  ub = fixef(fit_year_temp) + q * se
  )

# create both lb and ub model
fit_year_temp_lb <- fit_year_temp_ub <- fit_year_temp
fit_year_temp_lb@beta <- tab$lb
fit_year_temp_ub@beta <- tab$ub


# for each specie, make prediction and plot result
par(mfrow = c(3,3), las = 2, mar = c(4,3,2,0))
spec <- names(rmse)
for(i in 1:length(spec)){
  
  # predict upper and lower bound according to 
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
  
  # point estimate
  pred_pe <- (pred_ub - pred_lb)/2 + pred_lb
  
  
  # ylim range values
  ylim <- c(0, 
            max(
              Y$individualCount[Y$species == spec[i]],
              pred_lb, pred_ub
              )
  )
  
  plot(
    # Y$year[Y$species == s],
    X$year,
    c(
      Y$individualCount[Y$species == spec[i]],
      rep(NA, nrow(X) - sum(Y$species == spec[i]))
      ),
    type = 'l',
    ylim = ylim,
    xlab = '',
    ylab = paste(spec[i],'individual count'),
    main = paste0(i,': out of sample evolution'),
    xaxt = 'n'
  )
  
  
  # axis value
  if(i %% 9 > 6 | i %% 9 == 0) axis(1, at = X$year, labels = X$year, tick = FALSE)
  
  polygon(
    c(X$year, rev(X$year)),
    c(pred_lb, rev(pred_ub)),
    border = NA,
    col = ggplot2::alpha(2, .2)
  )
  
  lines(
    X$year,
    pred_pe,
    lwd = 1.5,
    col = 'red')
}























