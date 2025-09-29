#### Merge all data together ####

data2025<-data2025[,c(1:4,7,5,8:19,6,20,21)]
data2024 <- data2024[,-c(21:23)]
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

all.data <- all.data[,c(1:8,20,21,28,29,22:27,9:19)]

write.csv(all.data,file="./../clean_data/all_data_combined.csv",row.names = FALSE)