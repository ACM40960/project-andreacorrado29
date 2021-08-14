# SA Drought: consequences on species survival

## Abstract
Within this project, we explore and model the global temperature and rain over the last century with a system of ODE. Moreover, we keep track of different waterbirds species individual counts in South Africa over time. 
The interest of the project is to understand what the relation between the global climate change and the animal species count. We make use of the data available to train different models and validate them by specific diagnostic tools.
Once the different models validation, we want to make prediction up to a certain level of uncertainty. By making prediction for future years, we want to understand whether and when these species will extinguish.

## Overview 

What is happening to the global climate?

<br>
<img src="https://github.com/ACM40960/project-andreacorrado29/tree/master/figures/rain.gif" width=40%>
![Demo File](https://figures/rain.gif)

If you are interested in the full detailed study, check the report out `20205529_final.pdf` 

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

