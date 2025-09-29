#### fitting models with biomass as response

biomass <- data[!is.na(data$aboveground_mass),]

# biomass as a function of neighbors, microsite, and elevation
mass <- lm(log(aboveground_mass) ~ neighbor.number + patch + elevation,data=biomass)
summary(mass)
mass <- lm(log(aboveground_mass) ~ neighbor.biomass + patch + elevation,data=biomass)
summary(mass)

# biomass as a function of neighbors and environmental variables
biomass <- biomass[!is.na(biomass$mean.season.temp),]

all.env.biomass <- lm(log(aboveground_mass) ~ neighbor.biomass + mean.season.VWC + min.season.VWC + 
                        max.season.VWC + mean.season.temp + min.season.temp + max.season.temp, 
                      data=biomass)
mean.env.biomass <- lm(log(aboveground_mass) ~ neighbor.biomass + mean.season.VWC + mean.season.temp, 
                       data=biomass)
elevation.biomass <- lm(log(aboveground_mass) ~ neighbor.biomass + elevation,data=biomass)

# AIC to determine if elevation or environmental variables explain more variation
AIC(all.env.biomass,mean.env.biomass,elevation.biomass) # all env model is best

# how do neighborhood and environmental variables affect aboveground biomass?
ggplot(biomass,aes(x=neighbor.number,y=log(aboveground_mass),color=patch,shape=elevation))+
  geom_point()+
  labs(x="Number of neighbors",y="log(aboveground biomass)")
ggplot(biomass,aes(x=neighbor.biomass,y=log(aboveground_mass),color=patch,shape=elevation))+
  geom_point()+
  labs(x="Neighbor biomass",y="log(aboveground biomass)")
ggplot(biomass,aes(x=mean.season.VWC,y=log(aboveground_mass)))+
  geom_point()+
  labs(x="Mean soil moisture (%VWC)",y="log(aboveground biomass)")
ggplot(biomass,aes(x=mean.season.temp,y=log(aboveground_mass)))+
  geom_point()+
  labs(x="Mean soil temperature",y="log(aboveground biomass)")
ggplot(biomass,aes(x=patch,y=log(aboveground_mass)))+
  geom_boxplot()
ggplot(biomass,aes(x=elevation,y=log(aboveground_mass))) +
  geom_boxplot()
ggplot(biomass,aes(x=subplot,y=log(aboveground_mass))) +
  geom_boxplot()


#### FITNESS ####

fitness <- data[!is.na(data$fitness),]
fitness$fitness <- as.numeric(fitness$fitness)

# fitness as a function of neighbors, microsite, and elevation
fit <- lm(log(fitness) ~ neighbor.biomass + patch + elevation,data=fitness)
summary(fit)
fit <- lm(log(fitness) ~ neighbor.number + patch + elevation,data=fitness)
summary(fit)

# fitness as a function of neighbors and environmental variables
fitness.temp <- fitness[!is.na(fitness$mean.season.temp),]

all.env.fitness <- lm(log(fitness) ~ neighbor.biomass + mean.season.VWC + min.season.VWC + 
                        max.season.VWC + mean.season.temp + min.season.temp + max.season.temp, 
                      data=fitness.temp)
mean.env.fitness <- lm(log(fitness) ~ neighbor.biomass + mean.season.VWC + mean.season.temp, 
                       data=fitness.temp)
elevation.fitness <- lm(log(fitness) ~ neighbor.biomass + elevation,data=fitness.temp)

# AIC to determine if elevation or environmental variables explain more variation
AIC(all.env.fitness,mean.env.fitness,elevation.fitness) # elevation or all env variables model is best

ggplot(fitness,aes(x=neighbor.number,y=log(fitness),color=patch,shape=elevation))+
  geom_point()+
  labs(x="Number of neighbors",y="log(fitness)")
ggplot(fitness,aes(x=neighbor.biomass,y=log(fitness),color=patch,shape=elevation))+
  geom_point()+
  labs(x="Neighbor biomass",y="log(fitness)")
ggplot(fitness,aes(x=subplot,y=log(fitness)))+
  geom_boxplot()
ggplot(fitness,aes(x=elevation,y=log(fitness)))+
  geom_boxplot()
ggplot(fitness,aes(x=patch,y=log(fitness)))+
  geom_boxplot()

ggplot(fitness,aes(x=mean.season.VWC,y=log(fitness)))+
  geom_point()+
  labs(x="Mean soil moisture (%VWC)",y="log(fitness)")
ggplot(fitness,aes(x=mean.season.temp,y=log(fitness)))+
  geom_point()+
  labs(x="Mean soil temperature",y="log(fitness)")
