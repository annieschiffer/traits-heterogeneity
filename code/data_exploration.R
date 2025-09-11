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

glm(germination ~ neighbor.number + patch + elevation,data=data,family = "binomial")

# how do environmental variables affect germination?
env <- data[!is.na(data$mean.season.temp),]
all.env.variables <- glm(germination ~ neighbor.number + mean.season.VWC + min.season.VWC + max.season.VWC + mean.season.temp + min.season.temp + max.season.temp, 
    data=env,family="binomial")
mean.env.variables <- glm(germination ~ neighbor.number + mean.season.VWC + mean.season.temp, 
                         data=env,family="binomial")
just.elevation <- glm(germination ~ neighbor.number + elevation,data=env,family = "binomial")

# model selection on same data frame to see if elevation or the environmental variables explain more variance
AIC(all.env.variables,mean.env.variables,just.elevation) # all environmental variables explain more of the variance

#### BIOMASS ####

biomass <- data[!is.na(data$aboveground_mass),]
lm(log(aboveground_mass) ~ neighbor.number + patch + elevation,data=biomass)
lm(log(aboveground_mass) ~ neighbor.number + mean.season.VWC + min.season.VWC + max.season.VWC + mean.season.temp + min.season.temp + max.season.temp, 
    data=biomass)

# how do neighborhood and environmental variables affect aboveground biomass?
ggplot(biomass,aes(x=neighbor.number,y=log(aboveground_mass),color=patch,shape=elevation))+
  geom_point()
ggplot(biomass,aes(x=mean.season.VWC,y=log(aboveground_mass)))+
  geom_point()

  
#### TRAITS ####

## Phenology ##

# basic formatting
weekly.pheno <- weekly[-(which(is.na(weekly$germination))),c(1:8,10)]
weekly.pheno$phenophase[weekly.pheno$germination==1 & is.na(weekly.pheno$phenophase)]<-1
weekly.pheno<-na.omit(weekly.pheno)
weekly.pheno$date<-isoweek(as.Date(weekly.pheno$date,format = "%m/%d/%Y"))

# setting up to plot by subplot/competition treatment
pheno.treat <- weekly.pheno %>% group_by(date,subplot,phenophase) %>% summarize(abundance=n())
flower.treat <- pheno.treat[pheno.treat$phenophase==3,]
emerg.treat <- pheno.treat[pheno.treat$phenophase==1,]
senes.treat <- pheno.treat[pheno.treat$phenophase==4,]

# setting up to plot by patch
pheno.patch <- weekly.pheno %>% group_by(date,patch,phenophase) %>% summarize(abundance=n())
flower.patch <- pheno.patch[pheno.patch$phenophase==3,]
emerg.patch <- pheno.patch[pheno.patch$phenophase==1,]
senes.patch <- pheno.patch[pheno.patch$phenophase==4,]

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

# removing supplemental individuals
planted.mass<-above.mass[-(which(above.mass$subplot=="supp")),]
#planted.mass<-planted.mass[-grep("comp",planted.mass$shrub),]
fitness.mass <- merge(fitness,planted.mass, by=c("site","shrub","patch","subplot","species","id","elevation"))

# plotting
ggplot(above.mass,aes(x=species,y=aboveground_mass,fill=patch))+
  geom_boxplot()+
  labs(x="Species",y="Aboveground biomass (g)",fill="Patch")
ggplot(planted.mass,aes(x=species,y=aboveground_mass,fill=subplot))+
  geom_boxplot()+
  labs(x="Species",y="Aboveground biomass (g)",fill="Competition")
ggplot(above.mass,aes(x=species,y=aboveground_mass,fill=elevation))+
  geom_boxplot()+
  labs(x="Species",y="Aboveground biomass (g)",fill="Elevation")
ggplot(fitness.mass,aes(x=aboveground_mass,y=fitness,color=patch))+
  geom_point()

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

