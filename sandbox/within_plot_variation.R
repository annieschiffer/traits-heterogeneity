
# import data
data <- read.csv("./../clean_data/all_data_combined.csv")
data <- data[data$species=="BRTE",]

neighbors <- read.csv("./../clean_data/planted.natural.neighbors.csv")

#### differences between shrubs within a site? -----

t <- glm(germination ~ shrub,data=data,family="binomial")
summary(t)

t <- lm(log(aboveground_mass) ~ shrub,data=data)
summary(t)

ggplot(neighbors,aes(x=as.factor(shrub),y=planted.number))+
  geom_boxplot()
ggplot(neighbors,aes(x=as.factor(shrub),y=n.con.nb))+
  geom_boxplot()
# no differences in vital rates between shrubs within a site

#### differences between planted and natural neighbors? ------

neighbors$avg.mass.natural <- neighbors$mass.con.nb/neighbors$n.con.nb
neighbors$avg.mass.planted <- neighbors$planted.mass/neighbors$planted.number
nbavg<- rename(neighbors,"natural"="avg.mass.natural","planted"="avg.mass.planted")
nbavg <- pivot_longer(nbavg,cols = c("natural","planted"),names_to = "origin",
                       values_to = "average_mass")
nbavg <- nbavg[which(nbavg$average_mass >0),c(1:3,10,11)]

ggplot(nbavg,aes(x=as.factor(patch),y=log(average_mass),fill=origin))+
  geom_boxplot()
ggplot(nbavg,aes(x=as.factor(shrub),y=log(average_mass),fill=origin))+
  geom_boxplot()

t <- lm(log(average_mass) ~ shrub + patch,data=nbavg)
summary(t)
# natural cheatgrass always bigger than planted cheatgrass 
# but doesn't depend on within-site shrub variation or patch type



