## Make figures

## Fig 1
# multi-panel pictures of experimental setup

## Fig 2 
# fixed effects: germination rates & biomass (see drawing)


## Fig 3
# biomass reveals stress gradient hypothesis
# effect of conspecifics in low (favorable) vs. high (stressful) elevation

# all neighbors together
biomass$nb.number <- biomass$n.con.nb + biomass$n.het.nb
ggplot(biomass,aes(x=nb.number,y=log(aboveground_mass),color=elevation))+
  geom_point()+
  geom_smooth(method="lm")
# conspecific neighbors
ggplot(biomass,aes(x=n.con.nb,y=log(aboveground_mass),color=elevation))+
  geom_point()+
  geom_smooth(method="lm")
# heterospecific neighbors
ggplot(biomass,aes(x=n.het.nb,y=log(aboveground_mass),color=elevation))+
  geom_point()+
  geom_smooth(method="lm")

## Fig 4
# trait shifts in low vs. high and open vs. shrub (see drawing)




### OLD ---------

## trait distributions

ggplot(ltraits,aes(x=SLA,color=patch,linetype = subplot))+
  geom_freqpoly()
ggplot(rtraits,aes(x=SRL,color=patch,linetype = subplot))+
  geom_freqpoly()
ggplot(traits,aes(x=max.height,color=patch,linetype=subplot))+
  geom_freqpoly()

## trait boxplots
ggplot(ltraits,aes(x=patch,y=SLA,fill=subplot))+
  geom_boxplot()

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
