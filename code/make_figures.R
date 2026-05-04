## Make figures

#### load data and custom functions
source("fig_functions.R")
source("PCA.R")
load("./../outputs/2026-04-27/stan_fits/germ.fit.rda")
load("./../outputs/2026-04-27/stan_fits/SRL.fit.rda")
gdata <- read.csv("./../clean_data/all_data_combined.csv")
gdata <- gdata[gdata$species=="BRTE",]
pca_data <- read.csv("./../clean_data/pca_data.csv")
clim.cov <- read.csv("./../clean_data/fig1_climate_cover.csv")
SS.clim.data <- read.csv("./../clean_data/fig1_SS_climate.csv")

##### Fig 1: climate variables, SGH proof -----------

fig1 <- ggplot()+
  geom_point(data=clim.cov,aes(x=MAP,y=MAT.F,size = avg.cover),alpha=0.2)+
  geom_point(data=SS.clim.data, aes(x=MAP,y=MAT.F,color=elevation),size=4,shape=17)+
  labs(x="Mean Daily Precipitation (mm)",y="Mean Daily Temperature (F)",color="Site elevation",
       size="Cheatgrass % cover")+
  theme_minimal()+
  scale_color_manual(limits = c("high","low"),values=brewer.pal(11,"PuOr")[c(3,4)])+
  theme(axis.title = element_text(size=15),legend.title = element_text(size=15),axis.text = element_text(size=12),
        legend.text=element_text(size=12))

##### Fig 2: SGH with emergence --------------

gdata.sum <- gdata %>% group_by(elevation,patch,subplot) %>% summarize(prop.emerg = sum(germination)/n())
gdata.sum$elevation <- factor(gdata.sum$elevation,levels=c("low","high"))

fig2 <- ggplot(gdata.sum,aes(x=patch,y=prop.emerg,fill=as.factor(subplot)))+
  facet_grid(~elevation)+
  geom_bar(stat = "identity",position = position_dodge())+
  scale_fill_manual(limits=c("C","R"),labels=c("present","absent"),values = brewer.pal(11,"PuOr")[c(9,4)])+
  theme_minimal()+
  labs(y="Proportion established",fill="Competition",)+
  theme(axis.title.x = element_blank(),axis.title.y = element_text(size=15),axis.text.x = element_text(size=15),axis.text.y = element_text(size=12),
        legend.title = element_text(size=15),legend.text = element_text(size=12),strip.text = element_text(size=15))


#### Fig 3: PCA showing trait dispersion -----------------

fig3 <- ggplot() +
  geom_point(data=pca_data,aes(x=PC1,y=PC2,color=elevation,shape=patch),alpha=0.7,size=3)+
  theme_minimal()+
  geom_segment(data=loads,aes(x=0,y=0,xend=PC1*4,yend=PC2*4),color="black",
               arrow = arrow(length=unit(0.2,"cm")))+
  geom_text(data=loads,aes(x=PC1*4.2,y=PC2*4.5,label=c("SLA","SRL","height")),color="black",size=5)+
  scale_shape_manual(values=c(2,19))+
  scale_color_manual(values=brewer.pal(11,"PuOr")[c(3,10)])+
  labs(color="Elevation",shape="Patch")+
  theme(legend.title = element_text(size=15),legend.text = element_text(size=12),axis.title=element_text(size=15))


##### Fig 4: SRL across environments ----------------

# trait shifts in low vs. high and open vs. shrub - only SRL varied with anything

rtraits <- gdata[-which(is.na(gdata$length_cm)),]
rtraits <- rtraits[-which(is.na(rtraits$root_mass)),]
rtraits <- rtraits[-which(rtraits$root_mass==0),]
rtraits$SRL <- as.numeric(rtraits$length_cm)/as.numeric(rtraits$root_mass)

rsum <- rtraits[,c("patch","elevation","subplot","SRL")]
fig4 <- ggplot(rsum,aes(x=patch,y=log(SRL),fill=subplot)) +
  geom_boxplot()+
  scale_fill_manual(limits=c("C","R"),
                    labels=c("present","absent"),
                    values=brewer.pal(11,"PuOr")[c(2,8)])+
  labs(fill="Competition")+
  theme_classic()+
  theme(axis.title.x = element_blank(),legend.title = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=15),axis.text.y=element_text(size=12),legend.text = element_text(size=12))+
  annotate("text",label="long,\nthin roots",x=-0.1,y=11,size=5,color="#2D004B",fontface="bold")+
  annotate("text",label="short,\nthick roots",x=-0.1,y=8.5,size=5,color="#2D004B",fontface="bold")+
  coord_cartesian(clip="off",xlim=c(0.95,NA))+
  theme(plot.margin = margin(l=70))

# save output
if(!dir.exists(paste0("./../outputs/", Sys.Date(),"/"))) dir.create(paste0("./../outputs/", Sys.Date(),"/"))
ggsave(fig1,file=paste0("./../outputs/2026-04-27/fig1.jpeg"),height=5,width=7)
ggsave(fig2,file = paste0("./../outputs/2026-04-27/fig2.jpeg"),height = 6,width = 10)
ggsave(fig3,file = paste0("./../outputs/2026-04-27/fig3.jpeg"),height = 5,width = 7)
ggsave(fig4, file=paste0("./../outputs/2026-04-27/fig4.jpeg"),height=5,width=7)

##### Supplemental figures --------------

### Phenology figures

# formatting data
etraits <- all.data[-which(is.na(all.data$emergence)),]
etraits$emergence <- as.numeric(strftime(etraits$emergence, format = "%V"))

ftraits <- all.data[-which(is.na(all.data$flower)),]
ftraits$flower <- as.numeric(strftime(ftraits$flower, format = "%V"))

ecount <- etraits %>% group_by(patch,elevation,subplot,emergence) %>% summarize(emerged = n())
fcount <- ftraits %>% group_by(patch,elevation,subplot,flower) %>% summarise(flowered=n())

# plot phenology
eplot<- ggplot(ecount,aes(x=emergence,y=emerged,color=patch,linetype = subplot)) +
  facet_wrap(~elevation,axes="all",axis.labels="all")+
  geom_line()+
  theme_minimal()+
  labs(x="Week of emergence",y="Number emerged",linetype="Competitors",color="Patch")+
  theme(axis.title = element_text(size=15),
        legend.text = element_text(size=12),legend.title = element_text(size=15),
        strip.text = element_text(size=12))+
  scale_color_manual(values=brewer.pal(11,"PuOr")[c(4,9)])+
  scale_linetype_manual(limits=c("C","R"),labels=c("present","absent"),
                        values=c("solid","dashed"))
fplot <- ggplot(fcount,aes(x=flower,y=flowered,color=patch,linetype = subplot)) +
  geom_line()+
  facet_wrap(~elevation,axes="all",axis.labels="all")+
  theme_minimal()+
  labs(x="Week of flowering",y="Number flowered",linetype="Competitors",color="Patch")+
  theme(axis.title = element_text(size=15),
        legend.text = element_text(size=12),legend.title = element_text(size=15),
        strip.text = element_text(size=12))+
  scale_color_manual(values=brewer.pal(11,"PuOr")[c(4,9)])+
  scale_linetype_manual(limits=c("C","R"),labels=c("present","absent"),
                        values=c("solid","dashed"))

# save plot
ggsave(eplot,file = paste0("./../outputs/", Sys.Date(),"/supp_phenology.jpeg"),height = 6,width = 10)

### Soil conditions across elevation and patches

# calculate moisture differences
low.sm <- mean(all.data$mean.season.VWC[all.data$elevation=="low"]) 
high.sm <- mean(all.data$mean.season.VWC[all.data$elevation=="high"])
open.sm <- mean(all.data$mean.season.VWC[all.data$patch=="open"]) 
shrub.sm <- mean(all.data$mean.season.VWC[all.data$patch=="shrub"])

# calculate snowmelt date differences
spring.temps <- soil.temp[(soil.temp$date > "2024-03-01" & soil.temp$date < "2024-04-10") |
                         (soil.temp$date > "2025-03-01" & soil.temp$date < "2025-04-10"),]
snow.present <- spring.temps[which(spring.temps$temp < 35 & spring.temps$temp > 28),]
snow.absent <- spring.temps[-which(spring.temps$temp < 35 & spring.temps$temp > 28),]
snowmelt.dates <- snow.absent %>% group_by(year,elevation,patch) %>% summarize(snowmelt = min(date))
snowmelt.dates$snowmelt <- yday(snowmelt.dates$snowmelt)

low.date <- mean(snowmelt.dates$snowmelt[snowmelt.dates$elevation=="low"])
high.date <- mean(snowmelt.dates$snowmelt[snowmelt.dates$elevation=="high"])
open.date <- mean(snowmelt.dates$snowmelt[snowmelt.dates$patch=="open"])
shrub.date <- mean(snowmelt.dates$snowmelt[snowmelt.dates$patch=="shrub"])

# formatting soil moisture and temperature data to plot
moisture.data <- pivot_longer(data=all.data,cols = c("mean.season.VWC","min.season.VWC","max.season.VWC"),names_to = "anomoly",values_to = "percent_VWC")
temp.data <- pivot_longer(data=all.data,cols = c("mean.season.temp","min.season.temp","max.season.temp"),names_to = "anomoly",values_to = "degrees_C")

# plot soil moisture and temp
moisture.plot <- ggplot(moisture.data,aes(x=as.factor(elevation),y=percent_VWC,fill=as.factor(patch)))+
  facet_grid(~anomoly,labeller=labeller(anomoly=c("max.season.VWC"="maximum","mean.season.VWC"="mean",
                                                   "min.season.VWC"="minimum")))+
  geom_boxplot()+
  labs(x="Elevation",y="% VWC",fill="Patch")+
  theme_minimal()+
  theme(strip.text = element_text(size=12),axis.title = element_text(size=15),axis.text = element_text(size=12),
        legend.title = element_text(size=15),legend.text=element_text(size=12))
moisture.plot

