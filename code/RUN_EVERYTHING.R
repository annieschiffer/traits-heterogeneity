
# load packages
if (!require("tidyr")) install.packages("tidyr"); library(tidyr)
if (!require("dplyr")) install.packages("dplyr"); library(dplyr)
if (!require("rstudioapi")) install.packages("rstudioapi"); library(rstudioapi)
if (!require("ggplot2")) install.packages("ggplot2"); library(ggplot2)
if (!require("lubridate")) install.packages("lubridate"); library(lubridate)
if (!require("stringr")) install.packages("stringr"); library(stringr)
if (!require("RColorBrewer")) install.packages("RColorBrewer"); library(RColorBrewer)
if (!require("vegan")) install.packages("vegan"); library(vegan)
if (!require("bayesplot")) install.packages("bayesplot"); library(bayesplot)
if (!require("amadeus")) install.packages("amadeus"); library(amadeus)
if (!require("terra")) install.packages("terra"); library(terra)
if (!require("geodata")) install.packages("geodata"); library(geodata)
if (!require("sf")) install.packages("sf"); library(sf)
if (!require("rstan")) install.packages("rstan"); library(rstan)
if (!require("ggpubr")) install.packages("ggpubr"); library(ggpubr)

# set working directory
current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path)) # set working directory to location of this file

# set output directory
if(!dir.exists(paste0("./../outputs/",Sys.Date(),"/")))dir.create(paste0("./../outputs/",Sys.Date(),"/"))
if(!dir.exists(paste0("./../outputs/",Sys.Date(),"/stan_fits/")))dir.create(paste0("./../outputs/",Sys.Date(),"/stan_fits/"))

#### Data cleaning ####

source("clean_2024_data.R")

source("clean_2025_data.R")

source("clean_env_data.R")

source("clean_neighbor_data.R")

source("clean_germinate_data.R")

# merging to get clean data frame
source("merge_data.R")

#### Data analysis ####

## Aim 1: test stress gradient hypothesis

source("fit_germination.R")

source("fit_biomass.R")

## Aim 2: analyze trait shifts along gradient

source("fit_traits.R")

source("PCA.R")

#### Make figures ####

source("make_figures.R")

