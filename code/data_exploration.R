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

# import data
above.mass<-read.csv("./../raw_data/aboveground_mass_clean.csv")
root.SA<-read.csv("./../raw_data/root_SA_clean.csv")
root.mass<-read.csv("./../raw_data/root_mass_clean.csv")
supp.height<-read.csv("./../raw_data/supp_height_clean.csv")
seeds<-read.csv("./../raw_data/seeds_clean.csv")
weekly<-read.csv("./../raw_data/weekly.csv")
soil.moisture<-read.csv("./../raw_data/soil_moisture.csv")

#### GERMINATION AND FITNESS ####

# get fecundity by species and site/patch/treat
seeds<-seeds[-(grep("supp",seeds$subplot)),]

# get germination rates
weekly.germ<-weekly[-(which(is.na(weekly$germination))),2:8]
germ<-distinct(weekly.germ)
germ.remov<-germ[germ$subplot=="R",]
germ.rate<- germ.remov %>% group_by(site,patch,species) %>% summarize(germ.rate=sum(germination)/n())
germ.rate.sp<-germ %>% group_by(species) %>% summarize(germ.rate=sum(germination)/n())

# calculate fitness
fecundity<-seeds[seeds$seed_number>0,]
fitness<-left_join(fecundity,germ.rate.sp,by="species")
fitness<-fitness %>% mutate(fitness=seed_number*germ.rate)

# summarizing for plots
fitness.treat <- fitness %>% group_by(subplot) %>% summarize(fit.mean=mean(fitness),fit.se=(sd(fitness)/sqrt(n())))
fitness.patch <- fitness %>% group_by(patch) %>% summarize(fit.mean=mean(fitness),fit.se=(sd(fitness)/sqrt(n())))

# plots
ggplot(germ.rate[(germ.rate$site=="low_gate" | germ.rate$site=="low_north"),],
       aes(x=species,y=germ.rate,fill=patch))+
  geom_bar(stat = "identity",position = position_dodge())+
  theme_classic()+
  labs(x="Species",y="Germination rate",fill="Patch",title = "Germination in LOW sites")+
  scale_fill_brewer(palette="Dark2")+
  scale_y_continuous(limits=c(0,0.6),breaks=seq(0,0.6,by=0.1))
ggplot(germ.rate[(germ.rate$site=="high_gate" | germ.rate$site=="high_north" | germ.rate$site=="high_east"),],
       aes(x=species,y=germ.rate,fill=patch))+
  geom_bar(stat = "identity",position = position_dodge())+
  theme_classic()+
  labs(x="Species",y="Germination rate",fill="Patch",title = "Germination in HIGH sites")+
  scale_fill_brewer(palette="Dark2")+
  scale_y_continuous(limits=c(0,0.6),breaks=seq(0,0.6,by=0.1))
ggplot(fitness.treat,aes(x=subplot,y=fit.mean))+
  geom_pointrange(aes(ymin = (fit.mean-fit.se),ymax=(fit.mean+fit.se)))+
  labs(x="Competition treatment",y="Fitness")+
  theme_classic()+
  scale_x_discrete(limits=c("C","R"),labels=c("competition","removal"))
ggplot(fitness.patch,aes(x=patch,y=fit.mean))+
  geom_pointrange(aes(ymin = (fit.mean-fit.se),ymax=(fit.mean+fit.se)))+
  theme_classic()+
  labs(x="Patch",y="Fitness")

#### NEIGHBORHOOD ####

# reading in neighborhood data and formatting
nb <- read.csv("./../raw_data/neighborhood.csv")
nb$patch[nb$patch=="O"]<-"open"
nb$patch[nb$patch=="S"]<-"shrub"

# summarizing neighbor abundance and merging with fitness data
nb.abund <- nb %>% group_by(site,shrub,patch,species,id) %>% summarize(nb.abund=n())
nb.abund$subplot<-rep("C",nrow(nb.abund))
fitness.nb <- left_join(fitness,nb.abund)

# summarizing neighbor biomass and merging with fitness data
nb.mass <- above.mass[grep("comp",above.mass$shrub),]
nb.mass <- nb.mass %>% group_by(site,patch) %>% summarize(nb.mass=sum(aboveground_mass))
fitness.nb.mass <- merge(fitness[,-1],nb.mass, by=c("site","patch"))

# plotting
ggplot(na.omit(fitness.nb),aes(x=nb.abund,y=fitness,color=patch))+
  geom_point(size=4)+
  labs(x="Neighbor abundance",y="Fitness",color="Patch")+
  theme_classic()
ggplot(fitness.nb.mass,aes(x=nb.mass,y=fitness,color=patch))+
  geom_point(size=4)+
  labs(x="Neighbor biomass",y="Fitness",color="Patch")+
  theme_classic()
  
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

## Root traits

roots<- merge(root.mass[,-1],root.SA[,-1],by=c("site","shrub","patch","subplot","species","id"))
roots<- roots %>% mutate(SRL=sa_cm2/length_cm)

ggplot(roots,aes(x=species,y=SRL,fill=patch))+
  geom_boxplot()+
  labs(x="Species",y="Specific root length (cm2/g)",fill="Patch")
ggplot(roots,aes(x=species,y=SRL,fill=site))+  
  geom_boxplot()+
  labs(x="Species",y="Specific root length (cm2/g)",fill="Site")

## Aboveground biomass

planted.mass<-above.mass[-(which(above.mass$subplot=="supp")),]
#planted.mass<-planted.mass[-grep("comp",planted.mass$shrub),]

ggplot(planted.mass,aes(x=species,y=aboveground_mass,fill=patch))+
  geom_boxplot()
ggplot(planted.mass,aes(x=species,y=aboveground_mass,fill=subplot))+
  geom_boxplot()
ggplot(planted.mass,aes(x=species,y=aboveground_mass,fill=site))+
  geom_boxplot()

#### OLD ####

## Looking at all relationships
# 
# leaf.traits<- leaf.traits %>% 
#   mutate(date=NULL,quad=NULL,id=NULL)
# plot(leaf.traits)
# 
# seed.traits<- seed.traits %>%
#   mutate(date=NULL,quad=NULL,id=NULL)
# plot(seed.traits)
# 
# ## Looking at distributions of certain traits
# 
# hist(leaf.traits$SLA) # normally distributed
# hist(seed.traits$avg.seed.mass) # wonky but looks close to normal
# hist(complete.traits$height) # right skewed
# 
# # initial distribution fitting
# fit.SLA<-fitdist(leaf.traits$SLA,"norm")
# plot(fit.SLA)
# 
# fit.seedmass<-fitdist(seed.traits$avg.seed.mass,"norm")
# plot(fit.seedmass)
# 
# fit.height<-fitdist(complete.traits$height,"lnorm")
# plot(fit.height)
# 
# ## Plotting relationships between predictors and traits
# 
# # distance to shrub and traits
# plot(SLA ~ distance_shrub, data=leaf.traits, col=species)
# plot(avg.seed.mass ~ distance_shrub, col=species, data=seed.traits)
# plot(log(height) ~ distance_shrub, data=complete.traits)
# plot(phenophase ~ distance_shrub, data=complete.traits)
# 
# # aspect and traits
# boxplot(SLA ~ aspect, data=leaf.traits)
# boxplot(avg.seed.mass ~ aspect, data=seed.traits)
# boxplot(height ~ aspect,data=complete.traits)
# boxplot(phenophase ~ aspect, data=complete.traits)
# 
# # elevation and traits
# boxplot(SLA ~ elevation, data=leaf.traits)
# boxplot(avg.seed.mass ~ elevation, data=seed.traits)
# boxplot(height ~ elevation,data=complete.traits)
# boxplot(phenophase ~ elevation, data=complete.traits)
# 
# ## Fitting linear regressions
# 
# SLA.model<- lmer(SLA ~ distance_shrub*aspect + elevation + (1|species),data=leaf.traits)
# summary(SLA.model)
# 
# seedmass.model<-lmer(avg.seed.mass ~ distance_shrub*aspect + elevation + (1|species),data=seed.traits)
# summary(seedmass.model)
# 
# height.model<-lmer(log(height) ~ distance_shrub*aspect + elevation + (1|species),data=complete.traits)
# summary(height.model)

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
