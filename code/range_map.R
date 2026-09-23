### create figure/map to justify that sheep station is at cold & high end of cheatgrass's range

# read in Tyson's climate and cheatgrass cover dataset
range.data <- read.csv("./../raw_data/trainandtest11_10.csv")
range.data <- range.data[,c("longitude.x","latitude.x","pctcov","pctcovlog","pres2","year","Elevation","Aspect","Slope","annualtempmean")]
range.sum <- range.data %>% group_by(longitude.x,latitude.x,year) %>% summarize(avg.cov = mean(pctcov))

# get coordinates of Tyson's cheatgrass sites
coord <- unique(range.sum[,c("longitude.x","latitude.x")])
coord$id <- 1:nrow(coord)

# create a function to extract climate data from gridmet
getClimateData <- function(Year,Variable,VariableAbb,coord){
  
  # download data
  amadeus::download_data(
    dataset_name = "gridmet",
    variable = Variable,
    year = Year,
    directory_to_save = "./../raw_data/",
    acknowledgement = TRUE,
    download=TRUE
  )
  
  # convert data
  rast.data <- rast(paste0("./../raw_data/",VariableAbb,"/",VariableAbb,"_",Year,".nc"))
  # get precip data for the coordinates
  Coordinates <- vect(coord,geom=c("longitude.x","latitude.x"),crs="EPSG:4326")
  values <- extract(rast.data,Coordinates)
  #browser()
  # calculate annual mean precip and total annual precip
  out_data <- values %>%
    rowwise() %>%
    mutate(
      annual_mean = mean(c_across(-ID), na.rm = TRUE)
    )
  # format output data
  out_data <- out_data[,c("ID","annual_mean")]
  out_data$year <- Year
  out_data$variable <- Variable
  out_data <- bind_cols(coord,out_data)

  file.remove(paste0("./../raw_data/",VariableAbb,"/",VariableAbb,"_",Year,".nc"))
  return(out_data)
}

# get climate data for all years of Tyson's cheatgrass cover dataset
years <- c(2002:2004,2006:2016)
all.clim <- list()
for(i in 1:length(years)){
  mintemp <- getClimateData(years[i],"Minimum Near-Surface Air Temperature","tmmn",coord)
  maxtemp <- getClimateData(years[i],"Maximum Near-Surface Air Temperature","tmmx",coord)
  precip <- getClimateData(years[i],"Precipitation","pr",coord)
  clim.comb <- rbind(precip,mintemp,maxtemp)
  all.clim[[i]] <- pivot_wider(clim.comb,names_from = variable,values_from = annual_mean)
}

## get cliamte data for all years at Sheep Station

# create data frame of Sheep Station coordinates
SScoord <- data.frame(
  ID=1:6,
  longitude.x=c(-112.158785,-112.184183,-112.165726,-112.054741,-112.057581,-112.071352),
  latitude.x=c(44.227257,44.230490,44.236490,44.293420,44.302991,44.301233)
)

# extract gridmet climate data for the Sheep Station sites
SS.clim <- list()
for(i in 1:length(years)){
  precip <- getClimateData(years[i],"Precipitation","pr",SScoord)
  mintemp <- getClimateData(years[i],"Minimum Near-Surface Air Temperature","tmmn",SScoord)
  maxtemp <- getClimateData(years[i],"Maximum Near-Surface Air Temperature","tmmx",SScoord)
  #browser()
  clim.comb <- rbind(precip,mintemp,maxtemp)
  SS.clim[[i]] <- pivot_wider(clim.comb,names_from = variable,values_from = annual_mean)
}

# summarize and format for plotting - climate and cover data
all.clim.data <- do.call(rbind,all.clim)
clim.cov <- left_join(range.sum,all.clim.data,by=c("longitude.x","latitude.x","year"))
# rename columns
clim.cov <- clim.cov[,-c(5,6)]
names(clim.cov) <- c("lon","lat","year","avg.cover","MAP","MinAT","MaxAT")
clim.cov <- clim.cov %>% group_by(lon,lat) %>% summarize(MAP = mean(MAP),MinAT=mean(MinAT),MaxAT=mean(MaxAT),
                                                         MAT = mean((MinAT+MaxAT)/2),avg.cover=mean(avg.cover))
clim.cov$MAT.F <- (clim.cov$MAT - 273.15) * (9/5) + 32

# summarize and format for plotting - Sheep Station climate data
SS.clim.data <- do.call(rbind,SS.clim)
SS.clim.data <- SS.clim.data[,-c(1,4)]
names(SS.clim.data) <- c("lon","lat","year","MAP","MinAT","MaxAT")
SS.clim.data <- SS.clim.data %>% group_by(lon,lat) %>% summarize(MAP = mean(MAP),MinAT=mean(MinAT),
                                                                 MAT = mean((MinAT+MaxAT)/2),MaxAT=mean(MaxAT))
SS.clim.data$elevation <- c("low","low","low","high","high","high")
SS.clim.data$MAT.F <- (SS.clim.data$MAT - 273.15) * (9/5) + 32

# save output files
write.csv(clim.cov,"./../clean_data/fig1_climate_cover.csv",row.names=FALSE)
write.csv(SS.clim.data,"./../clean_data/fig1_SS_climate.csv",row.names=FALSE)

