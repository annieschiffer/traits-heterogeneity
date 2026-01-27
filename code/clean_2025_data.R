
#### 2025 data cleaning ####

# import and combine trait data
traits2025<-read.csv("./../raw_data/2025/traits.csv")
traits2025<- traits2025[1:114,]
root25 <- read.csv("./../raw_data/2025/root_sa25.csv")
traits2025 <- merge(traits2025,root25,all.x=TRUE,by.x="label",by.y="ID") # checked, no duplicates
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
phenology<- merge(emergence,flower,by=c("site","shrub","patch","subplot","species","id"),all.x=TRUE) # checked, no duplicates
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

# merging phenology and trait data - checked, no duplicates
traits2025 <- merge(phenology,traits2025,by=c("site","shrub","patch","subplot","species","id","elevation"),all.x = TRUE)
traits2025$year <- 2025

# quick data frame clean up
traits2025 <- traits2025[,c(1:11,15,18,14,16,17,12,13,21)]
names(traits2025)[12:13] <- c("aboveground_mass","length_cm")
traits2024 <- traits2024[,-c(14:15)]

rm(emergence,flower,fruit,phenology,root25,weekly.germ,weekly25)