# SA Drought: consequences on species survival

## Overview 
Interest to climate change has recently grown considerably. Here, we go through the main phenomenons such as rainfall amount and temperature to understand how their trend is evolving over time.

#### Rain & Temperature over time

<br>

<img src="https://github.com/ACM40960/project-andreacorrado29/blob/master/figures/rain.gif" width=40%>
<img src="https://github.com/ACM40960/project-andreacorrado29/blob/master/figures/temp.gif" width=40%>


#### Wildlife species survival 

Moreover, we are interested to the relationship this phenomenon has with the wildife survival. in particular, in *South Africa* there have been severe drought which have shown to considerably influence the animals death. 

<img src="https://github.com/ACM40960/project-andreacorrado29/blob/master/figures/count.gif" width=40%>

According to the current trend, what would happen in 50 years time?  check the report  `20205529_final.pdf` to find it out.

## Installation instruction

#### R
The whole project is coded in `R`. R is an open source language for statistical analysis.

For a Windows machine, R can be downloaded at the following link  <https://cran.r-project.org/bin/windows/base/>  and installed by running the executable `.exe`

On a Linux Debian based machine, R can be easily installed from the terminal by running
`sudo apt-get install r-base`

#### RStudio
*The RStudio IDE is a set of integrated tools designed to help you be more productive with R and Python. It includes a console, syntax-highlighting editor that supports direct code execution, and a variety of robust tools for plotting, viewing history, debugging and managing your workspace.*
You can download `RStudio` from <https://www.rstudio.com/products/rstudio/download/#download> and install it by running the executable.

#### R packages

In order to run the code you need to install different external packages. To do so, make sure your computer is connected to the web and run the script `pacakges_install.R` in your terminal:

- open the terminal

- get into the project folder `cd <project_folder_path>`
- run the script into the master folder `Rscript ./master/pacakges_install.R`

This script will install all the package need for the project that are not available yet in the local machine. If you wish to install other packages, you can do that by `install.packages('<package_name')` into the `R` console, where `<package_name>` is the name of the package you wish to install.


## Data Source

All the data analysed within this project, are available in this reporsitory at `./data/`. The occourrence data set has been split into 4 subsets to allows for git hub storage. The file `load_occourrence.R` allows to load the data combined into a unique object.
Those can be downloaded from the following web services, respectively waterbirds data and climate data:

- GBIF - Global Biodiversity Information Facility *Free and open access to biodiversity data* <www.gbif.org/occurrence/>
- World Bank Group *Climate Change Knowledge Portal* <https://climateknowledgeportal.worldbank.org/download-data>

## Result reproducibility

In order to reproduce the analysis and extend the results presented in this study, start `final_code.R` from the local folder

- `cd <project_folder_path>`
- `rstudio ./master/final_code.R`

The code is clearly commented, therefore, once you have opened the `final_code` into your environment, you can easily change the settings and play with the code according to your interest. Do not forget to contribute with the interesting findings you will have!

## Repository Files
Within this repository, you will have access to the data folder `data` which includes:
-  `occourrence<i>.txt` individual count data set been split to allows for git hub storage
-  `data_split_load.R`  the `R` script used to split the main data set into 4 subsets
-  `load_occourrence.R` `the R` script to load and combine all the `occourrence<i>.csv` into one data set
-  `pr_1901_2020_ZAF.csv`, raw precipitation data downloaded from the data source
-  `tas_1901_2020_ZAF.csv` raw temperature data downloaded from the data source
-  `Y_wide.csv`, individual count data by column result of the analysis
-  `climate_pred` climate data result of the analysis

The folder `master` contains:

- `20205529_final.pdf` the final report
- `20205529_final.Rmd` the `Rmarkdown` script used to produce the final report
- `sample.bib` the `biblatex` file containing the research references
- `final_code.R` the `R` source code for reproducibility
- `functions.R` the `R` script with useful function to be sourced
- `packages_install.R` the `R` script to install all the useful packages for the analysis

- `git_ACM40960.Rproj` the `RStudio` environment used for the project

The folder `figures` contains:

-  `*.gif` animation used within this documentation 

- `readme_fig.R`  R script to produce the `.png` and instruction to convert into `.gif`

  
