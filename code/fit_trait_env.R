### Fitting models with germination as response

# import data
data <- read.csv("./../clean_data/all_data_combined.csv")
data <- data[data$species=="BRTE",]

# only complete data set
traits <- data[,c("max.height","length_cm","root_mass","total_leaf_area","total_leaf_mass")]
complete <- data[which(complete.cases(traits)==TRUE),]

# calculate SLA
complete$total_leaf_area <- complete$total_leaf_area/100
complete$SLA <- as.numeric(complete$total_leaf_area)/as.numeric(complete$total_leaf_mass)

# calculate SRL
complete$SRL <- as.numeric(complete$length_cm)/as.numeric(complete$root_mass)
complete <- complete[-which(complete$SRL==Inf),]

# phenology
complete$emergence <- as.numeric(strftime(complete$emergence, format = "%V"))

# formatting site-year
complete$siteyear <- paste0(complete$site,complete$year)
complete$siteyear <- as.numeric(as.factor(complete$siteyear))
complete$env <- paste0(complete$elevation,"-",complete$patch)

# scaling
complete$max.height <- as.numeric(scale(complete$max.height))
complete$SLA <- as.numeric(scale(complete$SLA))
complete$SRL <- as.numeric(scale(complete$SRL))
complete$emergence <- as.numeric(scale(complete$emergence))

# biomass as a function of traits, environmental conditions, and their interaction
biomass.trait <- lmer(log(aboveground_mass) ~ env*SLA + env*SRL + env*max.height + env*emergence +
                (1|siteyear),data=complete)
summary(biomass.trait)
