### Plotting and testing at annual data

# run necessary scripts
source("data_cleaning.R")

# load packages
if (!require("tidyr")) install.packages("tidyr"); library(tidyr)
if (!require("dplyr")) install.packages("dplyr"); library(dplyr)
if (!require("lme4")) install.packages("lme4"); library(lme4)
if (!require("fitdistrplus")) install.packages("fitdistrplus"); library(fitdistrplus)

# looking at all relationships
leaf.traits<- leaf.traits %>% 
  mutate(date=NULL,quad=NULL,id=NULL)
plot(leaf.traits)

seed.traits<- seed.traits %>%
  mutate(date=NULL,quad=NULL,id=NULL)
plot(seed.traits)

# looking at distributions of certain traits
hist(leaf.traits$SLA) # normally distributed
hist(seed.traits$avg.seed.mass) # wonky but looks close to normal
hist(complete.traits$height) # right skewed

# looking at relationships between distance to shrub and certain traits
plot(SLA ~ distance_shrub, data=leaf.traits, col=species)
plot(avg.seed.mass ~ distance_shrub, col=species, data=seed.traits)
plot(height ~ distance_shrub, data=complete.traits)
plot(phenophase ~ distance_shrub, data=complete.traits)

# looking at relationships between aspect and certain traits
boxplot(SLA ~ aspect, data=leaf.traits)
boxplot(avg.seed.mass ~ aspect, data=seed.traits)
boxplot(height ~ aspect,data=complete.traits)
boxplot(phenophase ~ aspect, data=complete.traits)

# looking at relationships between elevation and certain traits
boxplot(SLA ~ elevation, data=leaf.traits)
boxplot(avg.seed.mass ~ elevation, data=seed.traits)
boxplot(height ~ elevation,data=complete.traits)
boxplot(phenophase ~ elevation, data=complete.traits)

# looking at relationships 

# preliminary linear models
summary(lm(SLA ~ distance_shrub,data=leaf.traits)) # marginally significant
summary(lm(avg.seed.mass ~ distance_shrub, data=seed.traits)) # not significant
summary(lm(height ~ distance_shrub, data=complete.traits)) # not significant
summary(lm(phenophase ~ distance_shrub, data=complete.traits)) # significant

summary(lm(SLA ~ aspect,data=leaf.traits)) # not significant
summary(lm(avg.seed.mass ~ aspect, data=seed.traits)) # not significant
summary(lm(height ~ aspect, data=complete.traits)) # not significant
summary(lm(phenophase ~ aspect, data=complete.traits)) # not significant

summary(lm(SLA ~ elevation,data=leaf.traits)) # very significant
summary(lm(avg.seed.mass ~ elevation, data=seed.traits)) # very significant
summary(lm(height ~ elevation, data=complete.traits)) # significant
summary(lm(phenophase ~ elevation, data=complete.traits)) # not significant
