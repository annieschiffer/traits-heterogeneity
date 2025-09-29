### fitting models with traits as response

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
