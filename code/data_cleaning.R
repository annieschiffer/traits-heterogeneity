### Cleaning annual trait data

# load packages
if (!require("tidyr")) install.packages("tidyr"); library(tidyr)
if (!require("dplyr")) install.packages("dplyr"); library(dplyr)
if (!require("rstudioapi")) install.packages("rstudioapi"); library(rstudioapi)
if (!require("ggplot2")) install.packages("ggplot2"); library(ggplot2)
if (!require("lubridate")) install.packages("lubridate"); library(lubridate)


# set working directory
current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path)) # set working directory to location of this file

# import data
traits23<-("./../raw_data/annual_traits")
above.mass<-("./../raw_data/aboveground_mass")
root.SA<-("./../raw_data/root_SA")
root.mass<-("./../raw_data/root_mass")
supp.height<-("./../raw_data/supp_height")
seeds<-("./../raw_data/seeds")
weekly<-("./../raw_data/weekly")
soil.moisture<-("./../raw_data/soil_moisture")

files<-c(above.mass,root.SA,root.mass,seeds)

for(ifile in 1:length(files)){
  d<-read.csv(paste0(files[ifile],".csv"))
  d<-separate(d,col="label",into=c("site","treat","id"),sep="-")
  d$treat<-gsub("supp","",d$treat)
  d<-separate(d,col="treat",into = c("shrub","subplot"),sep = "(?<=[0-9])(?=\\s?[A-Z])")
  d<-separate(d,col="subplot",into=c("patch","subplot"),sep="(?<=[A-Z])(?=\\s?[A-Z])")
  d$subplot[is.na(d$subplot)]<-"supp"
  d<-separate(d,col="id",into=c("species","id"),sep="(?<=[A-Z])(?=\\s?[0-9])")
  d$id[is.na(d$id)]<-1
  
  # renaming variables
  d$patch[d$patch=="O"]<-"open"
  d$patch[d$patch=="S"]<-"shrub"
  d$species[d$species=="AD"]<-"ALDE"
  d$species[d$species=="BT"]<-"BRTE"
  d$species[d$species=="CP"]<-"COPA"
  d$species[d$species=="PD"]<-"PODO"
  d$species[d$species=="suppAD"]<-"ALDE"
  d$species[d$species=="suppBT"]<-"BRTE"
  d$species[d$species=="suppCP"]<-"COPA"
  d$species[d$species=="suppPD"]<-"PODO"
  d$site[d$site=="HE"]<-"high_east"
  d$site[d$site=="HG"]<-"high_gate"
  d$site[d$site=="HN"]<-"high_north"
  d$site[d$site=="LG"]<-"low_gate"
  d$site[d$site=="LN"]<-"low_north"
  
  write.csv(d,paste0(files[ifile],"_clean.csv"))
}


#### OLD ####

# creating separate trait data frames
# leaf.traits<- traits %>% 
#   filter(!is.na(total_leaf_area)) %>% # remove any rows with NA in the total leaf area column
#   mutate(seed_number=NULL,total_seed_mass=NULL,root=NULL,notes=NULL,# remove irrelevant columns
#         site=as.factor(site), # making categorical variables factors
#         aspect=as.factor(aspect),
#         species=as.factor(species),
#         phenophase=as.factor(phenophase),
#         total_leaf_mass= as.numeric(total_leaf_mass), # make leaf columns numeric
#         total_leaf_area=as.numeric(total_leaf_area),
#         elevation=case_when(site=="low_north" ~ "low",
#                             site=="low_gate" ~ "low",
#                             site=="low_east" ~ "low",
#                             site=="high_bench" ~ "high",
#                             site=="high_10" ~ "high",
#                             site=="high_9" ~ "high"),
#         SLA=total_leaf_mass/total_leaf_area) # calculate specific leaf area in mg per cm2
# 
# seed.traits<- traits %>% 
#   filter(!is.na(total_seed_mass)) %>% # remove any rows with NA in the total seed mass column
#   mutate(total_leaf_mass=NULL,total_leaf_area=NULL,root=NULL,notes=NULL, # remove irrelevant columns
#          site=as.factor(site), # making categorical variables factors
#          aspect=as.factor(aspect),
#          species=as.factor(species),
#          phenophase=as.factor(phenophase),
#          total_seed_mass=as.numeric(total_seed_mass), # make seed columns numeric
#          seed_number=as.numeric(seed_number), 
#          elevation=case_when(site=="low_north" ~ "low",
#                              site=="low_gate" ~ "low",
#                              site=="low_east" ~ "low",
#                              site=="high_bench" ~ "high",
#                              site=="high_10" ~ "high",
#                              site=="high_9" ~ "high"),
#          avg.seed.mass=total_seed_mass/seed_number) # calculate average seed mass in mg
# 
# complete.traits <- traits %>% 
#   mutate(root=NULL,notes=NULL,
#          elevation=case_when(site=="low_north" ~ "low",
#                              site=="low_gate" ~ "low",
#                              site=="low_east" ~ "low",
#                              site=="high_bench" ~ "high",
#                              site=="high_10" ~ "high",
#                              site=="high_9" ~ "high")) %>%
#   drop_na(distance_shrub)
