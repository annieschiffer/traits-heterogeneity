### Plotting and testing at annual data

# load packages
if (!require("rstudioapi")) install.packages("rstudioapi"); library(rstudioapi)
if (!require("tidyr")) install.packages("tidyr"); library(tidyr)
if (!require("dplyr")) install.packages("dplyr"); library(dplyr)
if (!require("lme4")) install.packages("lme4"); library(lme4)
if (!require("fitdistrplus")) install.packages("fitdistrplus"); library(fitdistrplus)

# set working directory
current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path)) # set working directory to location of this file

# run necessary scripts
source("data_cleaning.R")

## Looking at all relationships

leaf.traits<- leaf.traits %>% 
  mutate(date=NULL,quad=NULL,id=NULL)
plot(leaf.traits)

seed.traits<- seed.traits %>%
  mutate(date=NULL,quad=NULL,id=NULL)
plot(seed.traits)

## Looking at distributions of certain traits

hist(leaf.traits$SLA) # normally distributed
hist(seed.traits$avg.seed.mass) # wonky but looks close to normal
hist(complete.traits$height) # right skewed

# initial distribution fitting
fit.SLA<-fitdist(leaf.traits$SLA,"norm")
plot(fit.SLA)

fit.seedmass<-fitdist(seed.traits$avg.seed.mass,"norm")
plot(fit.seedmass)

fit.height<-fitdist(complete.traits$height,"lnorm")
plot(fit.height)

## Plotting relationships between predictors and traits

# distance to shrub and traits
plot(SLA ~ distance_shrub, data=leaf.traits, col=species)
plot(avg.seed.mass ~ distance_shrub, col=species, data=seed.traits)
plot(log(height) ~ distance_shrub, data=complete.traits)
plot(phenophase ~ distance_shrub, data=complete.traits)

# aspect and traits
boxplot(SLA ~ aspect, data=leaf.traits)
boxplot(avg.seed.mass ~ aspect, data=seed.traits)
boxplot(height ~ aspect,data=complete.traits)
boxplot(phenophase ~ aspect, data=complete.traits)

# elevation and traits
boxplot(SLA ~ elevation, data=leaf.traits)
boxplot(avg.seed.mass ~ elevation, data=seed.traits)
boxplot(height ~ elevation,data=complete.traits)
boxplot(phenophase ~ elevation, data=complete.traits)

## Fitting linear regressions

SLA.model<- lmer(SLA ~ distance_shrub*aspect + elevation + (1|species),data=leaf.traits)
summary(SLA.model)

seedmass.model<-lmer(avg.seed.mass ~ distance_shrub*aspect + elevation + (1|species),data=seed.traits)
summary(seedmass.model)

height.model<-lmer(log(height) ~ distance_shrub*aspect + elevation + (1|species),data=complete.traits)
summary(height.model)

# # looking at relationships 
# 
# # preliminary linear models
# summary(lm(SLA ~ distance_shrub,data=leaf.traits)) # marginally significant
# summary(lm(avg.seed.mass ~ distance_shrub, data=seed.traits)) # not significant
# summary(lm(height ~ distance_shrub, data=complete.traits)) # not significant
# summary(lm(phenophase ~ distance_shrub, data=complete.traits)) # significant
# 
# summary(lm(SLA ~ aspect,data=leaf.traits)) # not significant
# summary(lm(avg.seed.mass ~ aspect, data=seed.traits)) # not significant
# summary(lm(height ~ aspect, data=complete.traits)) # not significant
# summary(lm(phenophase ~ aspect, data=complete.traits)) # not significant
# 
# summary(lm(SLA ~ elevation,data=leaf.traits)) # very significant
# summary(lm(avg.seed.mass ~ elevation, data=seed.traits)) # very significant
# summary(lm(height ~ elevation, data=complete.traits)) # significant
# summary(lm(phenophase ~ elevation, data=complete.traits)) # not significant
