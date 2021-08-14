# readme_figs.R

# script to create png figures

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

p + transition_reveal(year) 

# create plot for temperature
p <-
  ggplot(X, aes(x = year, y=rain_obs, color = group) ) + # add axis
  geom_line() + # add line 
  scale_color_manual(values=c("red")) + # add color 
  theme(legend.position="none") + # remove legend
  labs(x = "Year", y = 'Rain') # add labels
p

p + transition_reveal(year) 
