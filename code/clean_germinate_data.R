
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

# calculating the germination rates by species and adding to data frame
n.germ24 <- data2024 %>% group_by(species) %>% summarize(n.germ=sum(germination),
                                                         n.total=n())
n.germ24$germ.rate <- n.germ24$n.germ/n.germ24$n.total * 100

data2024 <- merge(data2024,n.germ24,by="species",all.x=TRUE)

# computing fitness based on aboveground mass and germination rate
data2024$fitness <- data2024$aboveground_mass * data2024$germ.rate

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

# calculating the germination rates by species and adding to data frame
n.germ25 <- data2025 %>% group_by(species) %>% summarize(n.germ=sum(germination),
                                                         n.total=n())
n.germ25$germ.rate <- n.germ25$n.germ/n.germ25$n.total * 100

germ.rate <- n.germ25$germ.rate

# computing fitness based on aboveground mass and germination rate
data2025$fitness <- data2025$aboveground_mass * germ.rate

rm(weekly25,tst)