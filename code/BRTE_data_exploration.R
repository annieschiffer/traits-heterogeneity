### Plotting and testing at annual data

# load packages
if (!require("rstudioapi")) install.packages("rstudioapi"); library(rstudioapi)
if (!require("tidyr")) install.packages("tidyr"); library(tidyr)
if (!require("dplyr")) install.packages("dplyr"); library(dplyr)
if (!require("lme4")) install.packages("lme4"); library(lme4)
if (!require("fitdistrplus")) install.packages("fitdistrplus"); library(fitdistrplus)
if (!require("ggpubr")) install.packages("ggpubr"); library(ggpubr)
if (!require("RColorBrewer")) install.packages("RColorBrewer"); library(RColorBrewer)

# set working directory
current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path)) # set working directory to location of this file

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


#### BIOMASS ####

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


#### TRAITS ####
traits <- data[,c("max.height","length_cm","root_mass","total_leaf_area","total_leaf_mass")]
traits <- data[which(complete.cases(traits)==TRUE),]

## leaf and root traits

# leaf traits
traits$total_leaf_area <- traits$total_leaf_area/100
traits$SLA <- as.numeric(traits$total_leaf_area)/as.numeric(traits$total_leaf_mass)

lfit <- lm(log(SLA) ~ neighbor.biomass + patch + elevation,data=traits)
summary(lfit)

ggplot(traits,aes(x=SLA,y=log(fitness),color=subplot,shape=patch))+
  geom_point()+
  scale_shape_manual(values = c(0,15))
ggplot(traits,aes(x=SLA,y=log(fitness),color=elevation))+
  geom_point()

# root traits
rtraits <- traits[-which(traits$root_mass==0),]
rtraits$SRL <- as.numeric(rtraits$length_cm)/as.numeric(rtraits$root_mass)

rfit <- lm(log(SRL) ~ neighbor.biomass + patch + elevation,data=rtraits)
summary(rfit)

ggplot(rtraits,aes(x=SRL,y=log(fitness),color=subplot,shape=patch))+
  geom_point()+
  scale_shape_manual(values = c(0,15))
ggplot(rtraits,aes(x=SRL,y=log(fitness),color=elevation))+
  geom_point()

## height

hfit <-lm(log(max.height) ~ neighbor.biomass + patch + elevation,data=traits)
summary(hfit)

ggplot(traits,aes(x=max.height,y=log(fitness),color=subplot,shape=patch))+
  geom_point()+
  scale_shape_manual(values = c(0,15))
ggplot(traits,aes(x=max.height,y=log(fitness),color=elevation))+
  geom_point()

## phenology

# emergence
etraits <- traits[!is.na(traits$emergence),]
temp<-as.Date(etraits$emergence, "%Y-%m-%d")
etraits$emergence<-format(temp, format="%m-%d")

ggplot(etraits,aes(x=emergence,y=log(fitness),color=as.factor(year)))+
  geom_point()

# flowering 
ptraits <- traits[!is.na(traits$flower),]
temp<-as.Date(ptraits$flower, "%Y-%m-%d")
ptraits$flower<-format(temp, format="%m-%d")

ggplot(ptraits,aes(x=flower,y=log(fitness),color=as.factor(year)))+
  geom_point()
ggplot(ptraits,aes(x=flower,color=patch))+
  geom_point()

## trait distributions

ggplot(traits,aes(x=SLA,color=patch,linetype = subplot))+
  geom_freqpoly()
ggplot(rtraits,aes(x=SRL,color=patch,linetype = subplot))+
  geom_freqpoly()
ggplot(traits,aes(x=max.height,color=patch,linetype=subplot))+
  geom_freqpoly()

# Elevation & environment

ggplot(data,aes(x=as.factor(year),y=mean.season.VWC,fill=as.factor(elevation)))+
  geom_boxplot()
ggplot(data,aes(x=as.factor(year),y=min.season.VWC,fill=as.factor(elevation)))+
  geom_boxplot()
ggplot(data,aes(x=as.factor(year),y=max.season.VWC,fill=as.factor(elevation)))+
  geom_boxplot()

ggplot(data,aes(x=as.factor(elevation),y=mean.season.temp))+
  geom_boxplot()
ggplot(data,aes(x=as.factor(elevation),y=min.season.temp))+
  geom_boxplot()
ggplot(data,aes(x=as.factor(elevation),y=max.season.temp))+
  geom_boxplot()
