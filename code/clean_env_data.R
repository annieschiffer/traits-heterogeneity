
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
