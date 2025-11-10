## Make figures

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
