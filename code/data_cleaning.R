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
above.mass24<-("./../raw_data/2024/aboveground_mass")
root.SA24<-("./../raw_data/2024/root_SA")
root.mass24<-("./../raw_data/2024/root_mass")
supp.height24<-("./../raw_data/2024/supp_height")
seeds24<-("./../raw_data/2024/seeds")
soil.moisture24<-("./../raw_data/2024/soil_moisture")
soil.moisture25<-("./../raw_data/2025/soil_moisture")

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
  
  d <- d %>% mutate(elevation=case_when(site=="low_gate" ~ "low",
                      site=="low_north" ~ "low",
                      site=="high_gate" ~ "high",
                      site=="high_east" ~ "high",
                      site=="high_north" ~ "high"))
  
  write.csv(d,paste0(files[ifile],"_clean.csv"))
}

supp.height<-read.csv("./../raw_data/supp_height.csv")
supp.height$patch[supp.height$patch=="O"]<-"open"
supp.height$patch[supp.height$patch=="S"]<-"shrub"
supp.height<-rename(supp.height,max.height=max_height..cm.)
write.csv(supp.height,"./../raw_data/supp_height_clean.csv")

# get germination rates
weekly24<-read.csv("./../raw_data/2024/weekly.csv")
weekly25<-read.csv("./../raw_data/2025/weekly.csv")

# germination rates for 2024
weekly.germ<-weekly24[-(which(is.na(weekly24$germination))),2:8]
n.intact<-nrow(distinct(weekly.germ[,1:6]))
germ <- distinct(weekly.germ[weekly.germ$germination==1,])
nrow(germ)/n.intact

# germination rates for 2025
weekly.germ25<-weekly25[-(which(is.na(weekly25$id))),1:6]
germ25<-nrow(distinct(weekly.germ25[,2:6]))
germ25/1152

# total germination rate across 2 years
(nrow(germ)+germ25)/(n.intact+1152)

## soil temp files

season<-c("summer24","winter23-24")
patch<-c("shrub","open","cage")
site<-c("LG","HN","LN","HG")
soil.temp<-data.frame()

for (iseason in 1:length(season)) {
  for(ipatch in 1:length(patch)){
    for(isite in 1:length(site)){
      if(file.exists(paste0("./../raw_data/",site[isite],patch[ipatch],"_",season[iseason],".csv"))){
      s.temp<-read.csv(paste0("./../raw_data/",site[isite],patch[ipatch],"_",season[iseason],".csv"))}else{next}
      s.temp<-s.temp[-1,-c(1,4:7)]
      colnames(s.temp)<-c("date","temp")
      s.temp <- separate(s.temp,col = "date",into = c("date","time"),sep=" ")
      
      s.temp$patch <- patch[ipatch]
      s.temp$site <- site[isite]
      s.temp$season <- season[iseason]
      
      soil.temp <- rbind(soil.temp,s.temp)
    }
  }
}

soil.temp$site[soil.temp$site=="HG"]<-"high_gate"
soil.temp$site[soil.temp$site=="HN"]<-"high_north"
soil.temp$site[soil.temp$site=="LG"]<-"low_gate"
soil.temp$site[soil.temp$site=="LN"]<-"low_north"

soil.temp <- soil.temp %>% mutate(elevation=case_when(site=="low_gate" ~ "low",
                                      site=="low_north" ~ "low",
                                      site=="high_gate" ~ "high",
                                      site=="high_east" ~ "high",
                                      site=="high_north" ~ "high"))
write.csv(soil.temp,file = "./../raw_data/soil_temperature_clean.csv")

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
