# readme_figs.R

# script to create png figures
# you can easily merge all the PNG into a unique GIF by running
# convert -delay 15 -loop 0 *.png <your_gif_name>.gif
  
rm(list = ls())
gc()
library(gganimate)

# load climate data
X <- read.csv('../data/climate_pred.csv')[,-1]
X <- X[complete.cases(X),]

# create fake group for coloring
group <- factor(rep('2', nrow(X)), levels = c(0,1,2))

# create plot for temperature
p <-
  ggplot(X, aes(x = year, y=temp_obs, color = group) ) + # add axis
  geom_line() + # add line 
  scale_color_manual(values=c("red")) + # add color 
  theme(legend.position="none") + # remove legend
  labs(x = "Year", y = 'Temperature') # add labels
p

# p + transition_reveal(year) 

# create plot for temperature
p <-
  ggplot(X, aes(x = year, y=rain_obs, color = group) ) + # add axis
  geom_line() + # add line 
  scale_color_manual(values=c("red")) + # add color 
  theme(legend.position="none") + # remove legend
  labs(x = "Year", y = 'Rain') # add labels
p

# p + transition_reveal(year) 


# y wide ----------------------------------------------------------------

Y_wide <- read.csv('../data/Y_wide.csv')[,-1]
P <- ncol(Y_wide)
X <- melt(Y_wide[,c(1,(P-1):P)], id.var = 'year')

#plot
p <-
  ggplot(X, aes(x=year,y=value,group=variable,colour=variable)) +
  geom_line(size = .5) +
  geom_line(aes(lty=variable)) +
  scale_y_log10(breaks=c(1,2,5,10,25)) +
  labs(x = "Year", y = 'Individual Count') # add labels
p

# p + transition_reveal(year)













