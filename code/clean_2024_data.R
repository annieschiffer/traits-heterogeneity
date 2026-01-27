### Cleaning and merging 2024 data ####

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
# removing any duplicates
#phenology <- phenology[-which(duplicated(phenology[,c(1:6)])),] # no duplicates

# merge everything together
aboveground.biomass <- traits.list[[1]]
root.length <- traits.list[[2]]
root.mass <- traits.list[[3]]
fecundity <- traits.list[[4]]
leaves<-traits.list[[5]]

# remove accidental duplicates
aboveground.biomass <- aboveground.biomass[-which(duplicated(aboveground.biomass[,c(1:6)])),]
root.length <- root.length[-which(duplicated(root.length[,c(1:6)])),]
#root.mass <- root.mass[-which(duplicated(root.mass[,c(1:6)]))] # none duplicated
fecundity <- fecundity[-which(duplicated(fecundity[,c(1:6)])),]
leaves <- leaves[-which(duplicated(leaves[,c(1:6)])),]

traits2024 <- merge(aboveground.biomass,root.length,all.x = TRUE,by=c("site","shrub","patch","subplot","species","id","elevation"))
traits2024 <- merge(traits2024,root.mass,all.x = TRUE,by=c("site","shrub","patch","subplot","species","id","elevation"))
traits2024 <- merge(traits2024,fecundity,all.x = TRUE,by=c("site","shrub","patch","subplot","species","id","elevation"))
traits2024 <- merge(traits2024,leaves[,c(1:8,32)],all.x=TRUE,by=c("site","shrub","patch","subplot","species","id","elevation"))

traits2024 <- traits2024[,-15]
traits2024 <- merge(phenology,traits2024,all.x=TRUE,by=c("site","shrub","patch","subplot","species","id","elevation"))
traits2024$year <- 2024

rm(d,weekly24,weekly.germ,phenology,traits.list,aboveground.biomass,root.length,root.mass,fecundity,leaves,emergence,flower,fruit)

