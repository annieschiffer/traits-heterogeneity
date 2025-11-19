### Fitting models with germination as response

# import data
data <- read.csv("./../clean_data/all_data_combined.csv")
data <- data[data$species=="BRTE",]
data$nb.number <- data$n.con.nb + data$n.het.nb
# formatting
data$siteyear <- paste0(data$site,data$year)
data$siteyear <- as.numeric(as.factor(data$siteyear))
# scaling
data$nb.number <- as.numeric(scale(data$nb.number))
data$max.season.VWC <- as.numeric(scale(data$max.season.VWC))
data$mean.season.VWC <- as.numeric(scale(data$mean.season.VWC))
data$min.season.VWC <- as.numeric(scale(data$min.season.VWC))


# germination as a function of neighbors, microsite, and elevation
germ <- glmer(germination ~ elevation*nb.number + patch*nb.number+ 
              (1|siteyear),data=data,family = "binomial")
summary(germ)

