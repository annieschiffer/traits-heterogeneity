#### fitting models with biomass as response

# import data
data <- read.csv("./../clean_data/all_data_combined.csv")
data <- data[data$species=="BRTE",]

# formatting
biomass <- data[!is.na(data$aboveground_mass),]
biomass$siteyear <- paste0(biomass$site,biomass$year)
biomass$siteyear <- as.numeric(as.factor(biomass$siteyear))
# scaling
biomass$n.con.nb <- as.numeric(scale(biomass$n.con.nb))
biomass$n.het.nb <- as.numeric(scale(biomass$n.het.nb))
biomass$mass.con.nb <- as.numeric(scale(biomass$mass.con.nb))
biomass$mass.het.nb <- as.numeric(scale(biomass$mass.het.nb))
biomass$max.season.VWC <- as.numeric(scale(biomass$max.season.VWC))
biomass$mean.season.VWC <- as.numeric(scale(biomass$mean.season.VWC))
biomass$min.season.VWC <- as.numeric(scale(biomass$min.season.VWC))

# model
mod <- lmer(log(aboveground_mass) ~ n.con.nb*patch + n.het.nb*patch + n.con.nb*elevation + n.het.nb*elevation + 
              (1|siteyear),data=biomass)
summary(mod)

# quick plots
boxplot(mean.season.VWC ~ patch,data = biomass)
boxplot(mean.season.VWC ~ elevation,data = biomass)
