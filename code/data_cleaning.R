### Cleaning annual trait data

# load packages
if (!require("tidyr")) install.packages("tidyr"); library(tidyr)
if (!require("dplyr")) install.packages("dplyr"); library(dplyr)
if (!require("rstudioapi")) install.packages("rstudioapi"); library(rstudioapi)
if (!require("ggplot2")) install.packages("ggplot2"); library(ggplot2)
if (!require("lubridate")) install.packages("lubridate"); library(lubridate)
if (!require("stringr")) install.packages("stringr"); library(stringr)

# set working directory
current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path)) # set working directory to location of this file

#### 2024 data cleaning ####

# import and clean up trait data
above.mass24<-("./../raw_data/2024/aboveground_mass")
root.SA24<-("./../raw_data/2024/root_SA")
root.mass24<-("./../raw_data/2024/root_mass")
seeds24<-("./../raw_data/2024/seeds")
leaves24<- ("./../raw_data/2024/leaf_traits")

files<-c(above.mass24,root.SA24,root.mass24,seeds24,leaves24)

traits.list <- list()
supp.list <- list()
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
  
  supp.list[[ifile]] <- d[which(d$subplot=="supp"),]
  traits.list[[ifile]] <- d[-which(d$subplot=="supp"),]
}

# import and clean up phenology data
weekly24<-read.csv("./../raw_data/2024/weekly.csv")
weekly.germ<-weekly24[-(which(is.na(weekly24$germination))),c(1:8,10,11)]
weekly.germ<-weekly.germ[which(weekly.germ$germination==1),]
weekly.germ$phenophase[weekly.germ$date < "4/24/2024"] <- 1
weekly.germ$height[is.na(weekly.germ$height)] <- 0

# creating data frames for each phenological event
emergence <- weekly.germ %>% group_by(site,shrub,patch,subplot,species,id) %>% summarise(emergence = min(date),
                                                                                         max.height = max(height))
flower <- weekly.germ[which(weekly.germ$phenophase > 2),]
flower <- flower %>% group_by(site,shrub,patch,subplot,species,id) %>% summarise(flower = min(date))
fruit <- weekly.germ[which(weekly.germ$phenophase > 3),]
fruit <- fruit %>% group_by(site,shrub,patch,subplot,species,id) %>% summarise(fruit = min(date))

# merging phenology data together
phenology<- merge(emergence,flower,by=c("site","shrub","patch","subplot","species","id"),all.x=TRUE)
phenology<- merge(phenology,fruit,by=c("site","shrub","patch","subplot","species","id"),all.x=TRUE)
phenology <- phenology %>% mutate(elevation=case_when(site=="low_gate" ~ "low",
                                      site=="low_north" ~ "low",
                                      site=="high_gate" ~ "high",
                                      site=="high_east" ~ "high",
                                      site=="high_north" ~ "high"))

# merge everything together
aboveground.biomass <- traits.list[[1]]
root.length <- traits.list[[2]]
root.mass <- traits.list[[3]]
fecundity <- traits.list[[4]]
leaves<-traits.list[[5]]

traits2024 <- merge(aboveground.biomass,root.length,all.x = TRUE,by=c("site","shrub","patch","subplot","species","id","elevation"))
traits2024 <- merge(traits2024,root.mass,all.x = TRUE,by=c("site","shrub","patch","subplot","species","id","elevation"))
traits2024 <- merge(traits2024,fecundity,all.x = TRUE,by=c("site","shrub","patch","subplot","species","id","elevation"))
traits2024 <- merge(traits2024,leaves[,c(1:8,32)],all.x=TRUE,by=c("site","shrub","patch","subplot","species","id","elevation"))

traits2024 <- traits2024[,-15]
traits2024 <- merge(phenology,traits2024,all.x=TRUE,by=c("site","shrub","patch","subplot","species","id","elevation"))

rm(weekly24,weekly.germ,phenology,traits.list,aboveground.biomass,root.length,root.mass,fecundity,leaves,emergence,flower,fruit)

#### 2025 data cleaning ####

# import and combine trait data
traits2025<-read.csv("./../raw_data/2025/traits.csv")
traits2025<- traits2025[1:114,]
root25 <- read.csv("./../raw_data/2025/root_sa25.csv")
traits2025 <- merge(traits2025,root25,all.x=TRUE,by.x="label",by.y="ID")
supp2025 <- traits2025[grep("supp",traits2025$label),]
traits2025 <- traits2025[-grep("supp",traits2025$label),]

# clean up rows and columns
traits2025 <- traits2025 %>% 
  mutate(parts = str_extract_all(traits2025$label, "[A-Za-z]+|[0-9]+")) %>% 
  unnest_wider(parts,names_sep = "_")
traits2025<-traits2025[,c(13:16,2:7,10:12)]
traits2025<-separate(traits2025,col="parts_3",into=c("patch","subplot"),sep="(?<=[A-Z])(?=\\s?[A-Z])")
names(traits2025)[1:5] <- c("site","shrub","patch","subplot","id")
traits2025$patch[traits2025$patch=="O"]<-"open"
traits2025$patch[traits2025$patch=="S"]<-"shrub"
traits2025$site[traits2025$site=="HE"]<-"high_east"
traits2025$site[traits2025$site=="HG"]<-"high_gate"
traits2025$site[traits2025$site=="HN"]<-"high_north"
traits2025$site[traits2025$site=="LG"]<-"low_gate"
traits2025$site[traits2025$site=="LN"]<-"low_north"
traits2025$site[traits2025$site=="LE"]<-"low_east"
traits2025 <- traits2025 %>% mutate(elevation=case_when(site=="low_gate" ~ "low",
                                      site=="low_north" ~ "low",
                                      site=="low_east" ~ "low",
                                      site=="high_gate" ~ "high",
                                      site=="high_east" ~ "high",
                                      site=="high_north" ~ "high"))
traits2025$species <- "BRTE"

# import and clean up phenology data
weekly25<-read.csv("./../raw_data/2025/weekly.csv")
weekly.germ<-weekly25[-(which(is.na(weekly25$phenophase))),c(1:6,8:9)]
weekly.germ$height[is.na(weekly.germ$height)] <- 0
weekly.germ$species <- "BRTE"

# creating data frames for each phenological event
emergence <- weekly.germ %>% group_by(site,shrub,patch,subplot,species,id) %>% summarise(emergence = min(date),
                                                                                         max.height = max(height))
flower <- weekly.germ[which(weekly.germ$phenophase > 2),]
flower <- flower %>% group_by(site,shrub,patch,subplot,species,id) %>% summarise(flower = min(date))
fruit <- weekly.germ[which(weekly.germ$phenophase > 3),]
fruit <- fruit %>% group_by(site,shrub,patch,subplot,species,id) %>% summarise(fruit = min(date))

# merging phenology data together
phenology<- merge(emergence,flower,by=c("site","shrub","patch","subplot","species","id"),all.x=TRUE)
phenology<- merge(phenology,fruit,by=c("site","shrub","patch","subplot","species","id"),all.x=TRUE)
phenology <- phenology %>% mutate(elevation=case_when(site=="low_gate" ~ "low",
                                                      site=="low_north" ~ "low",
                                                      site=="low_east" ~ "low",
                                                      site=="high_gate" ~ "high",
                                                      site=="high_east" ~ "high",
                                                      site=="high_north" ~ "high"))
# merging phenology and trait data








#### soil moisture and temperature data ####

## soil temp files

soil.moisture24<-("./../raw_data/2024/soil_moisture")

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

#### Germination data ####
weekly25<-read.csv("./../raw_data/2025/weekly.csv")

# germination rates for 2024
n.intact<-nrow(distinct(weekly.germ[,1:6]))
germ <- distinct(weekly.germ[weekly.germ$germination==1,])
nrow(germ)/n.intact

# germination rates for 2025
weekly.germ25<-weekly25[-(which(is.na(weekly25$id))),1:6]
germ25<-nrow(distinct(weekly.germ25[,2:6]))
germ25/1152

# total germination rate across 2 years
(nrow(germ)+germ25)/(n.intact+1152)

#### Supplemental individual traits ####

supp.height<-read.csv("./../raw_data/2024/supp_height.csv")
supp.height$patch[supp.height$patch=="O"]<-"open"
supp.height$patch[supp.height$patch=="S"]<-"shrub"
supp.height<-rename(supp.height,max.height=max_height..cm.)