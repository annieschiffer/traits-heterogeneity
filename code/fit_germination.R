### Fitting models with germination as response

# import data
data <- read.csv("./../clean_data/all_data_combined.csv")
data <- data[data$species=="BRTE",]

#### GERMINATION ####

# germination as a function of neighbors, microsite, and elevation
germ <- glm(germination ~ neighbor.number + patch + elevation,data=data,family = "binomial")
summary(germ)
germ <- glm(germination ~ neighbor.biomass + patch + elevation,data=data,family = "binomial")
summary(germ)

# germination as a function of neighbors and environment
env <- data[!is.na(data$mean.season.temp),]
all.env.variables <- glm(germination ~ neighbor.number + mean.season.VWC + min.season.VWC + 
                           max.season.VWC + mean.season.temp + min.season.temp + max.season.temp, 
                         data=env,family="binomial")
mean.env.variables <- glm(germination ~ neighbor.number + mean.season.VWC + mean.season.temp, 
                          data=env,family="binomial")
just.elevation <- glm(germination ~ neighbor.number + elevation,data=env,family = "binomial")

# model selection on same data frame to see if elevation or the environmental variables explain more variance
AIC(all.env.variables,mean.env.variables,just.elevation) # all environmental variables explain more of the variance
