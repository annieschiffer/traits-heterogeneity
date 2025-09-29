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
