### Cleaning annual trait data

# load packages
if (!require("tidyr")) install.packages("tidyr"); library(tidyr)
if (!require("dplyr")) install.packages("dplyr"); library(dplyr)
if (!require("rstudioapi")) install.packages("rstudioapi"); library(rstudioapi)

# set working directory
current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path)) # set working directory to location of this file

# import data
traits23<-read.csv("./../raw_data/annual_traits.csv")
above.mass<-read.csv("./../raw_data/aboveground_mass.csv")
root.SA<-read.csv("./../raw_data/root_SA.csv")
root.mass<-read.csv("./../raw_data/root_mass.csv")
supp.height<-read.csv("./../raw_data/supp_height.csv")
seeds<-read.csv("./../raw_data/seeds.csv")
weekly<-read.csv("./../raw_data/weekly.csv")
soil.moisture<-read.csv("./../raw_data/soil_moisture.csv")

# get germination rates
weekly.germ<-weekly[-(which(is.na(weekly$germination))),2:8]
germ<-distinct(weekly.germ)
germ.rate<- germ %>% group_by(site,patch,subplot) %>% summarize(germ.rate=sum(germination)/n())

# get fecundity by species and site/patch/treat
seeds<-seeds[-(grep("supp",seeds$Label)),]
seeds.sep1<-separate(seeds,col="Label",into=c("site","treat","id"),sep="-")
seeds.sep1$treat<-gsub("supp","",seeds.sep1$treat)
seeds.sep2<-separate(seeds.sep1,col="treat",into = c("shrub","subplot"),sep = "(?<=[0-9])(?=\\s?[A-Z])")
seeds.sep3<-separate(seeds.sep2,col="subplot",into=c("patch","subplot"),sep="(?<=[A-Z])(?=\\s?[A-Z])")
#seeds.sep3$subplot[is.na(seeds.sep3$subplot)]<-"supp"
seeds.sep3$subplot[is.na(seeds.sep3$subplot)]<-"C"
seeds.sep<-separate(seeds.sep3,col="id",into=c("species","id"),sep="(?<=[A-Z])(?=\\s?[0-9])")
seeds.sep$id[is.na(seeds.sep$id)]<-1

seeds.sep$patch[seeds.sep$patch=="O"]<-"open"
seeds.sep$patch[seeds.sep$patch=="S"]<-"shrub"
seeds.sep$species[seeds.sep$species=="AD"]<-"ALDE"
seeds.sep$species[seeds.sep$species=="BT"]<-"BRTE"
seeds.sep$species[seeds.sep$species=="CP"]<-"COPA"
seeds.sep$species[seeds.sep$species=="PD"]<-"PODO"
seeds.sep$site[seeds.sep$site=="HE"]<-"high_east"
seeds.sep$site[seeds.sep$site=="HG"]<-"high_gate"
seeds.sep$site[seeds.sep$site=="HN"]<-"high_north"
seeds.sep$site[seeds.sep$site=="LG"]<-"low_gate"
seeds.sep$site[seeds.sep$site=="LN"]<-"low_north"


fecundity<- seeds.sep %>% group_by(site,patch,subplot) %>% summarize(fecundity=sum(seed_number))
fitness<-left_join(germ.rate,fecundity,by=c("site","patch","subplot"))
fitness$fecundity[which(is.na(fitness$fecundity))]<-0
fitness <- fitness %>% mutate(fitness=fecundity*germ.rate)

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
  drop_na(distance_shrub)
