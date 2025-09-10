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
for(ifile in 1:length(files)){
  d<-read.csv(paste0(files[ifile],".csv"))
  d<-separate(d,col="label",into=c("site","treat","id"),sep="-")
  d$treat<-gsub("supp","",d$treat)
  d$treat<-gsub("comp","",d$treat)
  d<-separate(d,col="treat",into = c("shrub","subplot"),sep = "(?<=[0-9])(?=\\s?[A-Z])")
  d<-separate(d,col="subplot",into=c("patch","subplot"),sep="(?<=[A-Z])(?=\\s?[A-Z])")
  d <- d[-which(is.na(d$subplot)),]
  d<-separate(d,col="id",into=c("species","id"),sep="(?<=[A-Z])(?=\\s?[0-9])")
  
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
  
  traits.list[[ifile]] <- d
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
phenology$emergence<-as.Date(phenology$emergence,format="%m/%d/%Y")
phenology$flower<-as.Date(phenology$flower,format="%m/%d/%Y")
phenology$fruit<-as.Date(phenology$fruit,format="%m/%d/%Y")

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
traits2024$year <- 2024

rm(d,weekly24,weekly.germ,phenology,traits.list,aboveground.biomass,root.length,root.mass,fecundity,leaves,emergence,flower,fruit)

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
weekly.germ$patch[weekly.germ$patch=="O"]<-"open"
weekly.germ$patch[weekly.germ$patch=="S"]<-"shrub"

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
phenology$emergence<-as.Date(phenology$emergence,format="%m/%d/%y")
phenology$flower<-as.Date(phenology$flower,format="%m/%d/%y")
phenology$fruit<-as.Date(phenology$fruit,format="%m/%d/%y")

# merging phenology and trait data
traits2025 <- merge(phenology,traits2025,by=c("site","shrub","patch","subplot","species","id","elevation"),all.x = TRUE)
traits2025$year <- 2025

# quick data frame clean up
traits2025 <- traits2025[,c(1:11,15,18,14,16,17,12,13,21)]
names(traits2025)[12:13] <- c("aboveground_mass","length_cm")
traits2024 <- traits2024[,-c(14:15)]

rm(emergence,flower,fruit,phenology,root25,weekly.germ,weekly25)

#### soil moisture and temperature anomoly data ####

## clean up soil moisture data
soil.moisture24<-read.csv("./../raw_data/2024/soil_moisture.csv")
soil.moisture25 <- read.csv("./../raw_data/2025/soil_moisture.csv")
soil.moisture24$date <- as.Date(soil.moisture24$date,format = "%m/%d/%Y")
soil.moisture25$date <- as.Date(soil.moisture25$date,format = "%m/%d/%y")
# merge the two years together adn clean up names
soil.moisture <- na.omit(rbind(soil.moisture24,soil.moisture25))
soil.moisture$patch[soil.moisture$patch=="C"]<-"cage"
soil.moisture$patch[soil.moisture$patch=="O"]<-"open"
soil.moisture$patch[soil.moisture$patch=="S"]<-"shrub"
# get mean soil moisture for underneath shrubs vs. in open spaces
soil.moisture <- soil.moisture %>% group_by(date,site,shrub,patch) %>% summarize(mean.VWC = mean(VWC))
soil.moisture <- soil.moisture %>% mutate(year=case_when(date < "2025-01-01" ~ "2024",
                                                      date > "2025-01-01" ~ "2025"))

# calculate anamoly data
SM.anamolies <- soil.moisture %>% group_by(year,site,shrub,patch) %>% summarize(mean.season.VWC = mean(mean.VWC),
                                                                          min.season.VWC = min(mean.VWC),
                                                                          max.season.VWC = max(mean.VWC))
SM.anamolies<-SM.anamolies[-which(SM.anamolies$patch=="cage"),]


## clean up soil temperature data
season<-c("winter23-24","summer24","winter24-25","summer25")
patch<-c("shrub","open","cage")
site<-c("LG","HN","LN","HG")
soil.temp<-data.frame()

for (iseason in 1:length(season)) {
  for(ipatch in 1:length(patch)){
    for(isite in 1:length(site)){
      if(file.exists(paste0("./../raw_data/soil_temp/",site[isite],patch[ipatch],"_",season[iseason],".csv"))){
      s.temp<-read.csv(paste0("./../raw_data/soil_temp/",site[isite],patch[ipatch],"_",season[iseason],".csv"))}else{next}
      #browser()
      s.temp<-s.temp[-1,2:3]
      colnames(s.temp)<-c("date","temp")
      s.temp <- separate(s.temp,col = "date",into = c("date","time"),sep=" ")
      
      s.temp$patch <- patch[ipatch]
      s.temp$site <- site[isite]
      s.temp$season <- season[iseason]
      
      soil.temp <- rbind(soil.temp,s.temp)
    }
  }
}

# formatting
soil.temp$site[soil.temp$site=="HG"]<-"high_gate"
soil.temp$site[soil.temp$site=="HN"]<-"high_north"
soil.temp$site[soil.temp$site=="LG"]<-"low_gate"
soil.temp$site[soil.temp$site=="LN"]<-"low_north"
soil.temp <- soil.temp %>% mutate(elevation=case_when(site=="low_gate" ~ "low",
                                      site=="low_north" ~ "low",
                                      site=="high_gate" ~ "high",
                                      site=="high_east" ~ "high",
                                      site=="high_north" ~ "high"))
soil.temp$season[grep("winter",soil.temp$season)] <- "winter"
soil.temp$season[grep("summer",soil.temp$season)] <- "summer"
soil.temp$temp <- as.numeric(soil.temp$temp)
soil.temp$date <- as.Date(soil.temp$date,format = "%m/%d/%Y")
soil.temp <- soil.temp %>% mutate(year=case_when(date < "2024-01-01" ~ "2023",
                                                 date < "2025-01-01" & date >= "2024-01-01" ~ "2024",
                                                 date >= "2025-01-01" ~ "2025"))
# restricting dates within season
soil.temp <- soil.temp[(soil.temp$date > "2023-12-01" & soil.temp$date < "2024-04-10") |
                          (soil.temp$date > "2024-05-01" & soil.temp$date < "2024-07-02") |
                          (soil.temp$date > "2024-12-01" & soil.temp$date < "2025-04-10") |
                         (soil.temp$date > "2025-05-01" & soil.temp$date < "2025-06-26"),]

soil.temp$temp <- as.numeric(soil.temp$temp)


# calculate anamoly data
ST.anamolies <- soil.temp %>% group_by(year,season,site,patch) %>% summarize(mean.season.temp = mean(temp),
                                                                          min.season.temp = min(temp),
                                                                          max.season.temp = max(temp))

rm(soil.moisture24,soil.moisture25,soil.moisture,soil.temp,s.temp)

#### Neighbor data ####

# 2024 data

# counting number of neighbors
n24 <- read.csv("./../raw_data/2024/neighborhood.csv") # includes the ones that are planted
n24 <- n24 %>% group_by(site,shrub,patch,neighbor_sp) %>% summarize(n=n())
n24<- pivot_wider(n24,names_from = neighbor_sp,values_from = n)
n24 <- n24[-2,-5]
n24[is.na(n24)] <-0
n24$n.neighbors <- n24$COPA + n24$BRTE + n24$ALDE + n24$PODO
n24$patch[n24$patch=="O"]<-"open"
n24$patch[n24$patch=="S"]<-"shrub"

# getting natural neighbor biomass
nmass<- read.csv("./../raw_data/2024/aboveground_mass.csv") # DOESN'T include the planted ones
nmass <- nmass[grep("comp",nmass$label),]
nmass<-separate(nmass,col="label",into=c("site","treat","id"),sep="-")
nmass<-separate(nmass,col="treat",into = c("shrub","patch"),sep = "(?<=[0-9])(?=\\s?[A-Z])")
nmass<-separate(nmass,col="id",into=c("species","id"),sep="(?<=[A-Z])(?=\\s?[0-9])")
nmass$shrub <- gsub("comp","",nmass$shrub)
nmass$patch[nmass$patch=="O"]<-"open"
nmass$patch[nmass$patch=="S"]<-"shrub"
nmass$site[nmass$site=="HE"]<-"high_east"
nmass$site[nmass$site=="HG"]<-"high_gate"
nmass$site[nmass$site=="HN"]<-"high_north"
nmass$site[nmass$site=="LG"]<-"low_gate"
nmass$site[nmass$site=="LN"]<-"low_north"

# getting the planted neighbor biomass
nmass <- nmass %>% group_by(site,shrub,patch) %>% summarize(n.mass = sum(aboveground_mass))
planted.mass <- traits2024[-which(is.na(traits2024$aboveground_mass)),]
planted.mass <- planted.mass[which(planted.mass$subplot=="C"),]
planted.mass <- planted.mass %>% group_by(site,shrub,patch) %>% summarize(planted.mass = sum(aboveground_mass))

# merging the biomass of planted individuals and natural individuals
neighbors2024 <- merge(planted.mass,nmass,by=c("site","shrub","patch"),all.x=TRUE)
neighbors2024$n.mass[is.na(neighbors2024$n.mass)] <- 0
neighbors2024$neighbor.biomass <- neighbors2024$planted.mass + neighbors2024$n.mass

# cleaning up the merged data frame
neighbors2024 <- merge(neighbors2024,n24,by=c("site","shrub","patch"),all.x=TRUE)
neighbors2024 <- na.omit(neighbors2024)
neighbors2024 <- neighbors2024[-2,-c(4:5,7:10)]
names(neighbors2024)[5] <- "neighbor.number"
neighbors2024$year <- "2024"

rm(planted.mass,n24,nmass)

# 2025 data

# getting the number and biomass of natural neighbors
n25 <- read.csv("./../raw_data/2025/neighborhood.csv") # DOESN'T include the planted ones
n25$n.neighbors <- n25$PODO+n25$BRTE+n25$ALDE+n25$COPA
n25$n.biomass <- n25$PODO_mass + n25$BRTE_mass + n25$ALDE_mass + n25$COPA_mass
n25<-n25[,-c(1,5:13)]
n25$patch[n25$patch=="O"]<-"open"
n25$patch[n25$patch=="S"]<-"shrub"

# getting the number and biomass of planted neighbors
planted.mass <- traits2025[-which(is.na(traits2025$aboveground_mass)),]
planted.mass <- planted.mass[planted.mass$subplot=="C",]
planted.mass <- planted.mass %>% group_by(site,shrub,patch) %>% summarize(planted.mass = sum(aboveground_mass),
                                                                          planted.number = n())
neighbors2025 <- merge(n25,planted.mass,by=c("site","shrub","patch"),all.x=TRUE)
neighbors2025$planted.mass[is.na(neighbors2025$planted.mass)] <- 0
neighbors2025$planted.number[is.na(neighbors2025$planted.number)] <- 0

# remove planted individuals without neighbors
remove <- which(neighbors2025$planted.number==1 & neighbors2025$n.neighbors==0)
neighbors2025 <- neighbors2025[-remove,]

# add neighbors and planted ones together
neighbors2025$neighbor.biomass <- neighbors2025$planted.mass + neighbors2025$n.biomass
neighbors2025$neighbor.number <- neighbors2025$planted.number + neighbors2025$n.neighbors
neighbors2025 <- neighbors2025[,-c(4:7)]
neighbors2025$year <- "2025"

# bind 2024 and 2025 data together
neighbors <- rbind(neighbors2024,neighbors2025)

rm(planted.mass,n25,neighbors2024,neighbors2025)

#### Germination data ####

# marking germination for 2024 data
weekly24<-read.csv("./../raw_data/2024/weekly.csv")
weekly.germ<-weekly24[-(which(is.na(weekly24$germination))),c(1:8,10,11)]

# identify the 2024 plants that never germinated
never.germinated <- weekly.germ %>% group_by(site,shrub,patch,subplot,species,id) %>% summarize(germ=sum(germination))
never.germinated <- never.germinated[never.germinated$germ==0,]
never.germinated<-never.germinated[,-7]
never.germinated <- never.germinated %>% mutate(elevation=case_when(site=="low_gate" ~ "low",
                                          site=="low_north" ~ "low",
                                          site=="high_gate" ~ "high",
                                          site=="high_east" ~ "high",
                                          site=="high_north" ~ "high"))
never.germinated <- never.germinated %>% mutate(emergence = NA,
                                                max.height = NA,
                                                flower=NA,
                                                fruit=NA,
                                                aboveground_mass=NA,
                                                length_cm=NA,
                                                root_mass=NA,
                                                seed_number=NA,
                                                total_seed_mass=NA,
                                                total_leaf_area=NA,
                                                total_leaf_mass=NA,
                                                year="2024")
data2024 <- rbind(traits2024,never.germinated)
data2024 <- data2024 %>% mutate(germination=case_when(is.na(emergence)==TRUE ~ 0,
                                                      is.na(emergence)==FALSE ~ 1))

rm(weekly24,weekly.germ,never.germinated)

# marking germination for 2025
weekly25<-read.csv("./../raw_data/2025/weekly.csv")
tst <-expand.grid(site=c("low_gate","low_north","low_east","high_gate","high_north","high_east"),
            shrub=seq(1:4),
            patch=c("open","shrub"),
            subplot=c("C","R"),
            id=seq(1:12))
tst$year <- "2025"
tst$species <- "BRTE"
tst <- tst %>% mutate(elevation=case_when(site=="low_gate" ~ "low",
                                                      site=="low_north" ~ "low",
                                                      site=="low_east" ~ "low",
                                                      site=="high_gate" ~ "high",
                                                      site=="high_east" ~ "high",
                                                      site=="high_north" ~ "high"))

data2025 <- merge(tst,traits2025,by=c("site","shrub","patch","subplot","id","year","species","elevation"),all.x = TRUE)
data2025 <- data2025 %>% mutate(germination=case_when(is.na(emergence)==TRUE ~ 0,
                                          is.na(emergence)==FALSE ~ 1))

rm(weekly25,tst)

#### Merge all data together ####

data2025<-data2025[,c(1:4,7,5,8:19,6,20)]
all.data <- rbind(data2024,data2025)

# add in environmental variables
ST.anamolies <- ST.anamolies[ST.anamolies$season=="summer",]
ST.anamolies <- ST.anamolies[,-2]

all.data <- merge(all.data,SM.anamolies,by=c("year","site","shrub","patch"),all.x=TRUE)
all.data <- merge(all.data,ST.anamolies,by=c("year","site","patch"),all.x = TRUE)


# add in neighborhood data
neighbors$subplot <- "C"

all.data <- merge(all.data,neighbors,by=c("year","site","shrub","patch","subplot"),all.x = TRUE)
all.data$neighbor.number[is.na(all.data$neighbor.number)] <- 0
all.data$neighbor.biomass[is.na(all.data$neighbor.biomass)] <- 0

all.data <- all.data[,c(1:8,20,27,28,21:26,9:19)]

write.csv(all.data,file="./../clean_data/all_data_combined.csv")

#### Supplemental individual traits ####

# supp.height<-read.csv("./../raw_data/2024/supp_height.csv")
# supp.height$patch[supp.height$patch=="O"]<-"open"
# supp.height$patch[supp.height$patch=="S"]<-"shrub"
# supp.height<-rename(supp.height,max.height=max_height..cm.)