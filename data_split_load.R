# deal with data splitting

# --------------------------------------------------------------------------- #
#                                data source                                  #
# --------------------------------------------------------------------------- #

# citation
# GBIF.org (09 June 2021) GBIF Occurrence Download https://doi.org/10.15468/dl.zy5hhv
link_data_query <- 'https://www.gbif.org/occurrence/'

# load in full data
path_to_full_file <- '~/Dropbox/ucd_MATH_PROJ/DROUGHT/occurrence.txt' 
dat <- read.delim(path_to_full_file)

# split in four parts
N <- nrow(dat)
N1 <- 1 : (N/4)
N2 <- (floor(N/4) + 1) :  (N/2)
N3 <- (N/2 + 1) : (floor(N/4) * 3)
N4 <- (floor(N/4) * 3 + 1) : N

dat1 <- dat[N1, ]
dat2 <- dat[N2, ]
dat3 <- dat[N3, ]
dat4 <- dat[N4, ]

# check file existence and write files
files_available <- list.files('./data/')
if (!('occurrence1.txt' %in% files_available)) write.csv(dat1, './data/occurrence1.txt')
if (!('occurrence2.txt' %in% files_available)) write.csv(dat1, './data/occurrence2.txt')
if (!('occurrence3.txt' %in% files_available)) write.csv(dat2, './data/occurrence3.txt')
if (!('occurrence4.txt' %in% files_available)) write.csv(dat2, './data/occurrence4.txt')
