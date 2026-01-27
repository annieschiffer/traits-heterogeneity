#### Merge all data together ####

data2025<-data2025[,c(1:4,7,5,8:19,6,20,21)]
data2024 <- data2024[,-c(21:23)]
all.data <- rbind(data2024,data2025) # no duplicates

# add in environmental variables
ST.anamolies <- ST.anamolies[ST.anamolies$season=="summer",]
ST.anamolies <- ST.anamolies[,-2]

all.data <- merge(all.data,SM.anamolies,by=c("year","site","shrub","patch"),all.x=TRUE)
all.data <- merge(all.data,ST.anamolies,by=c("year","site","patch"),all.x = TRUE)

# add in neighborhood data
neighbors$subplot <- "C"

all.data <- merge(all.data,neighbors,by=c("year","site","shrub","patch","subplot"),all.x = TRUE)
all.data$mass.het.nb[is.na(all.data$mass.het.nb)] <- 0
all.data$n.het.nb[is.na(all.data$n.het.nb)] <- 0
all.data$mass.con.nb[is.na(all.data$mass.con.nb)] <- 0
all.data$n.con.nb[is.na(all.data$n.con.nb)] <- 0

# quick conversion of dates to day of year
# all.data$emergence <- yday(all.data$emergence)
# all.data$flower <- yday(all.data$flower)
# all.data$fruit <- yday(all.data$fruit)
all.data <- all.data[,c(1:8,20,28:31,22:27,9:19)] # no duplicates

write.csv(all.data,file="./../clean_data/all_data_combined.csv",row.names = FALSE)