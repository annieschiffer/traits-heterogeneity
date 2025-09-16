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

#### GERMINATION ####

# germination as a function of neighbors, microsite, and elevation
germ <- glmer(germination ~ neighbor.number + patch + elevation + (1 | species),data=data,family = "binomial")
summary(germ)

# germination as a function of neighbors and environment
env <- data[!is.na(data$mean.season.temp),]
all.env.variables <- glmer(germination ~ neighbor.number + mean.season.VWC + min.season.VWC + 
                             max.season.VWC + mean.season.temp + min.season.temp + max.season.temp + (1 | species), 
    data=env,family="binomial")
mean.env.variables <- glmer(germination ~ neighbor.number + mean.season.VWC + mean.season.temp + (1 | species), 
                         data=env,family="binomial")
just.elevation <- glmer(germination ~ neighbor.number + elevation + (1 | species),data=env,family = "binomial")

# model selection on same data frame to see if elevation or the environmental variables explain more variance
AIC(all.env.variables,mean.env.variables,just.elevation) # all environmental variables explain more of the variance

#### BIOMASS ####

biomass <- data[!is.na(data$aboveground_mass),]

# biomass as a function of neighbors, microsite, and elevation
mass <- lmer(log(aboveground_mass) ~ neighbor.number + patch + elevation + (1 | species),data=biomass)
summary(mass)

# biomass as a function of neighbors and environmental variables
biomass <- biomass[!is.na(biomass$mean.season.temp),]

all.env.biomass <- lmer(log(aboveground_mass) ~ neighbor.number + mean.season.VWC + min.season.VWC + 
                        max.season.VWC + mean.season.temp + min.season.temp + max.season.temp + (1 | species), 
    data=biomass)
mean.env.biomass <- lmer(log(aboveground_mass) ~ neighbor.number + mean.season.VWC + mean.season.temp + (1 | species), 
                       data=biomass)
elevation.biomass <- lmer(log(aboveground_mass) ~ neighbor.number + elevation + (1 | species),data=biomass)

# AIC to determine if elevation or environmental variables explain more variation
AIC(all.env.biomass,mean.env.biomass,elevation.biomass) # elevation model is best

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


#### FITNESS ####

fitness <- data[!is.na(data$fitness),]
fitness$fitness <- as.numeric(fitness$fitness)

# biomass as a function of neighbors, microsite, and elevation
fit <- lmer(log(fitness) ~ neighbor.biomass + patch + elevation + (1 | species),data=fitness)
summary(fit)

ggplot(fitness,aes(x=neighbor.biomass,y=log(fitness)))+
  geom_point()



  
#### TRAITS ####
traits <- data[,c("max.height","length_cm","root_mass","total_leaf_area","total_leaf_mass")]
traits <- data[which(complete.cases(traits)==TRUE),]

## Phenology ##

# plotting subplot/competition treatment
ggplot(emerg.treat,aes(x=date,y=abundance,color=as.factor(subplot)))+
  geom_point()+
  geom_line()+
  labs(x="Week",y="Number in vegetation phenophase",color="Competition")
ggplot(flower.treat,aes(x=date,y=abundance,color=as.factor(subplot)))+
  geom_point()+
  geom_line()+
  labs(x="Week",y="Number in flowering phenophase",color="Competition")
ggplot(senes.treat,aes(x=date,y=abundance,color=as.factor(subplot)))+
  geom_point()+
  geom_line()+
  labs(x="Week",y="Number in senescence phenophase",color="Competition")

# plotting patch
ggplot(emerg.patch,aes(x=date,y=abundance,color=as.factor(patch)))+
  geom_point()+
  geom_line()+
  labs(x="Week",y="Number in vegetation phenophase",color="Patch")
ggplot(flower.patch,aes(x=date,y=abundance,color=as.factor(patch)))+
  geom_point()+
  geom_line()+
  labs(x="Week",y="Number in flowering phenophase",color="Patch")
ggplot(senes.patch,aes(x=date,y=abundance,color=as.factor(patch)))+
  geom_point()+
  geom_line()+
  labs(x="Week",y="Number in senescence phenophase",color="Patch")

## Root traits ##

# combining data frames and formatting
roots<- merge(root.mass[,-1],root.SA[,-1],by=c("site","shrub","patch","subplot","species","id","elevation"))
roots<- roots %>% mutate(SRL=sa_cm2/length_cm)
fitness.root <- merge(fitness,roots, by=c("site","shrub","patch","subplot","species","id","elevation"))

# plotting
ggplot(roots,aes(x=species,y=SRL,fill=patch))+
  geom_boxplot()+
  labs(x="Species",y="Specific root length (cm2/g)",fill="Patch")
ggplot(roots,aes(x=species,y=SRL,fill=elevation))+
  geom_boxplot()+
  labs(x="Species",y="Specific root length (cm2/g)",fill="Elevation")
ggplot(roots,aes(x=species,y=SRL,fill=subplot))+
  geom_boxplot()+
  labs(x="Species",y="Specific root length (cm2/g)",fill="Competition")
ggplot(fitness.root,aes(x=SRL,y=fitness,color=patch))+
  geom_point()

## Aboveground biomass ##

## Max height ##

# extracting height data from weekly data and combining with supps
weekly.height<-weekly[-(which(is.na(weekly$height))),c(2:7,11)]
height<-distinct(weekly.height)
height <- height %>% group_by(site,shrub,patch,subplot,species,id) %>% summarize(max.height=max(height))
all.height<- rbind(height[,-6],supp.height[,-1])
all.height <- all.height %>% mutate(elevation=case_when(site=="low_gate" ~ "low",
                                                site=="low_north" ~ "low",
                                                site=="high_gate" ~ "high",
                                                site=="high_east" ~ "high",
                                                site=="high_north" ~ "high"))
fitness.height <- merge(fitness,height, by=c("site","shrub","patch","subplot","species","id"))

# plotting
ggplot(all.height,aes(x=species,y=max.height,fill=patch))+
  geom_boxplot()+
  labs(x="Species",y="Max height (cm)",fill="Patch")
ggplot(height,aes(x=species,y=max.height,fill=subplot))+
  geom_boxplot()+
  labs(x="Species",y="Max height (cm)",fill="Competition")
ggplot(all.height,aes(x=species,y=max.height,fill=elevation))+
  geom_boxplot()+
  labs(x="Species",y="Max height (cm)",fill="Elevation")
ggplot(fitness.height,aes(x=max.height,y=fitness,color=patch))+
  geom_point()


# plotting soil temps
soil.temp.avg<-soil.temp.avg[!(soil.temp.avg$patch=="cage"),]

ggplot(soil.temp.avg[soil.temp.avg$date > "2024-01-01" & soil.temp.avg$date < "2024-07-10",],aes(x=date,y=avg.temp,color=patch))+
  #geom_point(size=3)+
  geom_line(aes(linetype=elevation))+
  theme_classic()+
  scale_color_brewer(palette="Dark2")+
  labs(x="Date",y="Temperature (ºC)",color="Patch",linetype="Elevation",title="Soil temperature over the season")

