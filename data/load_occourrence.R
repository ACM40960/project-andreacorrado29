
## load occourrence subsets and combine into a unique object
dat0 <- rbind(
  read.csv("occurrence1.txt", row.names=1),
  read.csv("occurrence2.txt", row.names=1),
  read.csv("occurrence3.txt", row.names=1),
  read.csv("occurrence4.txt", row.names=1)
)
