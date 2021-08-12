# SA Drought: consequences on species survival

## Abstract
Within this project, we explore and model the global temperature and rain over the last century with a system of ODE. Moreover, we keep track of different waterbirds species individual counts in South Africa over time. 
The interest of the project is to understand what the relation between the global climate change and the animal species count. We make use of the data available to train different models and validate them by specific diagnostic tools.
Once the different models validation, we want to make prediction up to a certain level of uncertainty. By making prediction for future years, we want to understand whether and when these species will extinguish.

## Installation instruction
The whole project is coded in `R`. R is an open source language for statistical analysis.

For a Windows machine, R can be downloaded at the following link  <https://cran.r-project.org/bin/windows/base/>  and installed by running the executable `.exe`

On a Linux Debian based machine, R can be easily installed from the terminal by running
`sudo apt-get install r-base`

In order to run the code you need to install different external packages. 
In order to do so, run the script `pacakges_install.R` script in your terminal

- get into the project folder `cd <project_folder_path>`
- run the script into the master folder `Rscript ./master/pacakges_install.R`

This script will install all the package need for the project that are not available yet in the local machine. If you wish to install other packages, you can do that by `install.packages('<package_name')` into the `R` console, where `<package_name>` is the name of the package you wish to install.


## Data Source

All the data analysed within this project, are available in this reporsitory at `./data/`.

Those can be downloaded from the following web services, respectively waterbirds data and climate data:

- GBIF - Global Biodiversity Information Facility *Free and open access to biodiversity data* <www.gbif.org/occurrence/>
- World Bank Group *Climate Change Knowledge Portal* <https://climateknowledgeportal.worldbank.org/download-data>

## Result reproducibility

