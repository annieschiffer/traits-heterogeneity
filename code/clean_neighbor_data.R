#### Neighbor data ####

# 2024 data

# counting number of neighbors
nb <- read.csv("./../raw_data/2024/neighborhood.csv") # includes the ones that are planted
nb <- nb %>% group_by(site,shrub,patch,neighbor_sp) %>% summarize(n=n())
nb <- nb[-2,]
nb$patch[nb$patch=="O"]<-"open"
nb$patch[nb$patch=="S"]<-"shrub"
nb$site[6] <- "high_north"

# getting neighbor biomass
all.mass<- read.csv("./../raw_data/2024/aboveground_mass.csv") # DOES include the planted ones
all.mass <- all.mass[-c(grep("supp",all.mass$label),grep("R",all.mass$label)),]
all.mass<-separate(all.mass,col="label",into=c("site","treat","id"),sep="-")
all.mass<-separate(all.mass,col="treat",into = c("shrub","patch"),sep = "(?<=[0-9])(?=\\s?[A-Z])")
all.mass<-separate(all.mass,col="id",into=c("species","id"),sep="(?<=[A-Z])(?=\\s?[0-9])")
all.mass$shrub <- gsub("comp","",all.mass$shrub)
all.mass$patch[all.mass$patch=="SC"] <- "S"
all.mass$patch[all.mass$patch=="OC"] <- "O"
all.mass$patch[all.mass$patch=="O"]<-"open"
all.mass$patch[all.mass$patch=="S"]<-"shrub"
all.mass$site[all.mass$site=="HE"]<-"high_east"
all.mass$site[all.mass$site=="HG"]<-"high_gate"
all.mass$site[all.mass$site=="HN"]<-"high_north"
all.mass$site[all.mass$site=="LG"]<-"low_gate"
all.mass$site[all.mass$site=="LN"]<-"low_north"
all.mass$species[all.mass$species=="CP"] <- "COPA"
all.mass$species[all.mass$species=="BT"] <- "BRTE"
all.mass$species[all.mass$species=="AD"] <- "ALDE"
all.mass$species[all.mass$species=="PD"] <- "PODO"

# get mass and count for each species
all.mass <- all.mass %>% group_by(site,shrub,patch,species) %>% summarize(n.mass = sum(aboveground_mass),
                                                                          n = n())
# remove the ones where BRTE is alone (no neighbors)
remove <- all.mass %>% group_by(site,shrub,patch) %>% summarize(n.emerged=sum(n),
                                                                    species=unique(species))
remove <- remove[which(remove$species=="BRTE" & remove$n.emerged==1),]
all.mass <- all.mass[-c(16,17),]

# merge number of neighbors from mass dataset and neighborhood dataset
mass.num <- all.mass[,-5]
names(nb) <- c("site","shrub","patch","species","n")
nb$shrub <- as.character(nb$shrub)
add.n.nb <- rbind(nb,mass.num)
nb.count <- add.n.nb %>% group_by(site,shrub,patch,species) %>% summarize(n.nb=max(n))

# pivot wider
nb.count<- pivot_wider(nb.count,names_from = species,values_from = n.nb)
nb.count[is.na(nb.count)] <-0
nb.count$n.het.nb <- nb.count$COPA + nb.count$ALDE + nb.count$PODO
nb.count$n.con.nb <- nb.count$BRTE
nb.count <- nb.count[,-c(4:7)]

nb.mass <- pivot_wider(all.mass,names_from = species,values_from = n.mass)
nb.mass[is.na(nb.mass)] <- 0
nb.mass$mass.het.nb <- nb.mass$COPA + nb.mass$ALDE + nb.mass$PODO
nb.mass$mass.con.nb <- nb.mass$BRTE
nb.mass <- nb.mass[,-c(4:8)]

# merge mass and number
nb2024 <- merge(nb.mass,nb.count,by=c("site","shrub","patch"),all=TRUE)
nb2024 <- nb2024[-14,]
nb2024$year <- "2024"

rm(nb,all.mass,mass.num,remove,add.n.nb,nb.mass)

# 2025 data

# getting the number and biomass of natural neighbors
n25 <- read.csv("./../raw_data/2025/neighborhood.csv") # DOESN'T include the planted ones
n25$n.het.nb <- n25$PODO+n25$ALDE+n25$COPA
n25$n.con.nb <- n25$BRTE
n25$mass.het.nb <- n25$PODO_mass + n25$ALDE_mass + n25$COPA_mass
n25$mass.con.nb <- n25$BRTE_mass
n25<-n25[,-c(1,5:13)]
n25$patch[n25$patch=="O"]<-"open"
n25$patch[n25$patch=="S"]<-"shrub"

# getting the number and biomass of planted neighbors
planted.mass <- traits2025[-which(is.na(traits2025$aboveground_mass)),]
planted.mass <- planted.mass[planted.mass$subplot=="C",]
planted.mass <- planted.mass %>% group_by(site,shrub,patch) %>% summarize(planted.mass = sum(aboveground_mass),
                                                                          planted.number = n())
nb2025 <- merge(n25,planted.mass,by=c("site","shrub","patch"),all.x=TRUE)
nb2025$planted.mass[is.na(nb2025$planted.mass)] <- 0
nb2025$planted.number[is.na(nb2025$planted.number)] <- 0

# remove planted individuals without neighbors
remove <- which(nb2025$planted.number==1 & nb2025$n.con.nb==0 &nb2025$n.het.nb==0)
nb2025 <- nb2025[-remove,]

# add neighbors and planted ones together
nb2025$mass.con.nb <- nb2025$planted.mass + nb2025$mass.con.nb
nb2025$n.con.nb <- nb2025$planted.number + nb2025$n.con.nb
nb2025 <- nb2025[,-c(8,9)]
nb2025 <- nb2025[,c(1,2,3,6,7,4,5)]
nb2025$year <- "2025"

# bind 2024 and 2025 data together
neighbors <- rbind(nb2024,nb2025)

rm(planted.mass,n25,nb2024,nb2025)
