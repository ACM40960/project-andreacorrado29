# deal with data splitting

# --------------------------------------------------------------------------- #
#                                data source                                  #
# --------------------------------------------------------------------------- #

# citation
# GBIF.org (09 June 2021) GBIF Occurrence Download https://doi.org/10.15468/dl.zy5hhv
link <- 'https://www.gbif.org/occurrence/'

# load in full data
path_to_full_file <- '~/Dropbox/ucd_MATH_PROJ/DROUGHT/occurrence.txt' 
dat <- read.delim(path_to_full_file)

# split in halves
N <- nrow(dat)
dat1 <- dat[ 1 : (N/2) , ]
dat2 <- dat[( N/2 + 1):N, ]

# check file existence and write files
files_available <- list.files('./data/')
if (!('occurrence1.txt' %in% files_available)) write.csv(dat1, './data/occurrence1.txt')
if (!('occurrence2.txt' %in% files_available)) write.csv(dat2, './data/occurrence2.txt')

