### Cleaning annual trait data

# load packages
if (!require("tidyr")) install.packages("tidyr"); library(tidyr)
if (!require("dplyr")) install.packages("dplyr"); library(dplyr)
if (!require("rstudioapi")) install.packages("rstudioapi"); library(rstudioapi)
if (!require("ggplot2")) install.packages("ggplot2"); library(ggplot2)
if (!require("lubridate")) install.packages("lubridate"); library(lubridate)


# set working directory
current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path)) # set working directory to location of this file

# import data
traits23<-read.csv("./../raw_data/annual_traits.csv")
above.mass<-read.csv("./../raw_data/aboveground_mass.csv")
root.SA<-read.csv("./../raw_data/root_SA.csv")
root.mass<-read.csv("./../raw_data/root_mass.csv")
supp.height<-read.csv("./../raw_data/supp_height.csv")
seeds<-read.csv("./../raw_data/seeds.csv")
weekly<-read.csv("./../raw_data/weekly.csv")
soil.moisture<-read.csv("./../raw_data/soil_moisture.csv")

#### GERMINATION AND FITNESS ####

# get fecundity by species and site/patch/treat
seeds<-seeds[-(grep("supp",seeds$Label)),]
seeds.sep1<-separate(seeds,col="Label",into=c("site","treat","id"),sep="-")
seeds.sep1$treat<-gsub("supp","",seeds.sep1$treat)
seeds.sep2<-separate(seeds.sep1,col="treat",into = c("shrub","subplot"),sep = "(?<=[0-9])(?=\\s?[A-Z])")
seeds.sep3<-separate(seeds.sep2,col="subplot",into=c("patch","subplot"),sep="(?<=[A-Z])(?=\\s?[A-Z])")
#seeds.sep3$subplot[is.na(seeds.sep3$subplot)]<-"supp"
seeds.sep3$subplot[is.na(seeds.sep3$subplot)]<-"C"
seeds.sep<-separate(seeds.sep3,col="id",into=c("species","id"),sep="(?<=[A-Z])(?=\\s?[0-9])")
seeds.sep$id[is.na(seeds.sep$id)]<-1

# renaming variables
seeds.sep$patch[seeds.sep$patch=="O"]<-"open"
seeds.sep$patch[seeds.sep$patch=="S"]<-"shrub"
seeds.sep$species[seeds.sep$species=="AD"]<-"ALDE"
seeds.sep$species[seeds.sep$species=="BT"]<-"BRTE"
seeds.sep$species[seeds.sep$species=="CP"]<-"COPA"
seeds.sep$species[seeds.sep$species=="PD"]<-"PODO"
seeds.sep$site[seeds.sep$site=="HE"]<-"high_east"
seeds.sep$site[seeds.sep$site=="HG"]<-"high_gate"
seeds.sep$site[seeds.sep$site=="HN"]<-"high_north"
seeds.sep$site[seeds.sep$site=="LG"]<-"low_gate"
seeds.sep$site[seeds.sep$site=="LN"]<-"low_north"

# get germination rates
weekly.germ<-weekly[-(which(is.na(weekly$germination))),2:8]
germ<-distinct(weekly.germ)
germ.rate<- germ %>% group_by(patch) %>% summarize(germ.rate=sum(germination)/n())

# calculate fitness
fecundity<-seeds.sep[seeds.sep$seed_number>0,]
fitness<-left_join(fecundity,germ.rate,by="species")
fitness<-fitness %>% mutate(fitness=seed_number*germ.rate)

# summarizing for plots
fitness.treat <- fitness %>% group_by(subplot) %>% summarize(fit.mean=mean(fitness),fit.se=(sd(fitness)/sqrt(n())))
fitness.patch <- fitness %>% group_by(patch) %>% summarize(fit.mean=mean(fitness),fit.se=(sd(fitness)/sqrt(n())))

# plots
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

# summarizing and merging with fitness data
nb.sum <- nb %>% group_by(site,shrub,patch,species,id) %>% summarize(nb.abund=n())
nb.sum$shrub<-as.character(nb.sum$shrub)
nb.sum$id<-as.character(nb.sum$id)
nb.sum$subplot<-rep("C",nrow(nb.sum))
fitness.nb <- left_join(fitness,nb.sum)

# plotting
ggplot(na.omit(fitness.nb),aes(x=nb.abund,y=fitness,color=patch))+
  geom_point(size=4)+
  labs(x="Neighbor abundance",y="Fitness",color="Patch")+
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



#### OLD ####

# creating separate trait data frames
leaf.traits<- traits %>% 
  filter(!is.na(total_leaf_area)) %>% # remove any rows with NA in the total leaf area column
  mutate(seed_number=NULL,total_seed_mass=NULL,root=NULL,notes=NULL,# remove irrelevant columns
        site=as.factor(site), # making categorical variables factors
        aspect=as.factor(aspect),
        species=as.factor(species),
        phenophase=as.factor(phenophase),
        total_leaf_mass= as.numeric(total_leaf_mass), # make leaf columns numeric
        total_leaf_area=as.numeric(total_leaf_area),
        elevation=case_when(site=="low_north" ~ "low",
                            site=="low_gate" ~ "low",
                            site=="low_east" ~ "low",
                            site=="high_bench" ~ "high",
                            site=="high_10" ~ "high",
                            site=="high_9" ~ "high"),
        SLA=total_leaf_mass/total_leaf_area) # calculate specific leaf area in mg per cm2

seed.traits<- traits %>% 
  filter(!is.na(total_seed_mass)) %>% # remove any rows with NA in the total seed mass column
  mutate(total_leaf_mass=NULL,total_leaf_area=NULL,root=NULL,notes=NULL, # remove irrelevant columns
         site=as.factor(site), # making categorical variables factors
         aspect=as.factor(aspect),
         species=as.factor(species),
         phenophase=as.factor(phenophase),
         total_seed_mass=as.numeric(total_seed_mass), # make seed columns numeric
         seed_number=as.numeric(seed_number), 
         elevation=case_when(site=="low_north" ~ "low",
                             site=="low_gate" ~ "low",
                             site=="low_east" ~ "low",
                             site=="high_bench" ~ "high",
                             site=="high_10" ~ "high",
                             site=="high_9" ~ "high"),
         avg.seed.mass=total_seed_mass/seed_number) # calculate average seed mass in mg

complete.traits <- traits %>% 
  mutate(root=NULL,notes=NULL,
         elevation=case_when(site=="low_north" ~ "low",
                             site=="low_gate" ~ "low",
                             site=="low_east" ~ "low",
                             site=="high_bench" ~ "high",
                             site=="high_10" ~ "high",
                             site=="high_9" ~ "high")) %>%
  drop_na(distance_shrub)
