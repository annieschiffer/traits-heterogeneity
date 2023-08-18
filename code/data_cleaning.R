### Cleaning annual trait data

# load packages
if (!require("tidyr")) install.packages("tidyr"); library(tidyr)
if (!require("dplyr")) install.packages("dplyr"); library(dplyr)
if (!require("rstudioapi")) install.packages("rstudioapi"); library(rstudioapi)

# set working directory
current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path)) # set working directory to location of this file

# import data
traits<-read.csv("./../raw_data/annual_traits.csv")

# creating separate trait data frames
leaf.traits<- traits %>% 
  filter(!is.na(total_leaf_area)) %>% # remove any rows with NA in the total leaf area column
  mutate(seed_number=NULL,total_seed_mass=NULL,root=NULL,notes=NULL,# remove irrelevant columns
        site=as.factor(site), # making categorical variables factors
        aspect=as.factor(aspect),
        species=as.factor(species),
        phenophase=as.factor(phenophase),
        total_leaf_mass= as.numeric(total_leaf_mass), # make leaf columns numeric
        total_leaf_area=as.numeric(total_leaf_area),
        elevation=case_when(site=="low_north" ~ "low",
                            site=="low_gate" ~ "low",
                            site=="low_east" ~ "low",
                            site=="high_bench" ~ "high",
                            site=="high_10" ~ "high",
                            site=="high_9" ~ "high"),
        SLA=total_leaf_mass/total_leaf_area) # calculate specific leaf area in mg per cm2

seed.traits<- traits %>% 
  filter(!is.na(total_seed_mass)) %>% # remove any rows with NA in the total seed mass column
  mutate(total_leaf_mass=NULL,total_leaf_area=NULL,root=NULL,notes=NULL, # remove irrelevant columns
         site=as.factor(site), # making categorical variables factors
         aspect=as.factor(aspect),
         species=as.factor(species),
         phenophase=as.factor(phenophase),
         total_seed_mass=as.numeric(total_seed_mass), # make seed columns numeric
         seed_number=as.numeric(seed_number), 
         elevation=case_when(site=="low_north" ~ "low",
                             site=="low_gate" ~ "low",
                             site=="low_east" ~ "low",
                             site=="high_bench" ~ "high",
                             site=="high_10" ~ "high",
                             site=="high_9" ~ "high"),
         avg.seed.mass=total_seed_mass/seed_number) # calculate average seed mass in mg

complete.traits <- traits %>% 
  mutate(root=NULL,notes=NULL,
         elevation=case_when(site=="low_north" ~ "low",
                             site=="low_gate" ~ "low",
                             site=="low_east" ~ "low",
                             site=="high_bench" ~ "high",
                             site=="high_10" ~ "high",
                             site=="high_9" ~ "high")) %>%
  na.omit()
