## Make figures

#### load data and custom functions
source("fig_functions.R")
source("PCA.R")
# whole dataset
gdata <- read.csv("./../clean_data/all_data_combined.csv")
gdata <- gdata[gdata$species=="BRTE",]
# formatting
gdata$siteyear <- paste0(gdata$year,gdata$site,gdata$shrub)
gdata$siteyear <- as.numeric(as.factor(gdata$siteyear))
bdata <- gdata[!is.na(gdata$aboveground_mass),]
# other datasets
pca_data <- read.csv("./../clean_data/pca_data.csv")
clim.cov <- read.csv("./../clean_data/fig1_climate_cover.csv")
SS.clim.data <- read.csv("./../clean_data/fig1_SS_climate.csv")

##### Fig 1: climate variables, SGH proof -----------

fig1 <- ggplot()+
  geom_point(data=clim.cov,aes(x=MAP,y=MAT.F,size = avg.cover),alpha=0.2)+
  geom_point(data=SS.clim.data, aes(x=MAP,y=MAT.F,fill=elevation),size=6,color="black",shape=23)+
  labs(x="Mean Daily Precipitation (mm)",y="Mean Daily Temperature (F)",fill="Site\nelevation",
       size="Cheatgrass\n% cover")+
  theme_minimal()+
  scale_fill_manual(limits = c("high","low"),values=brewer.pal(11,"PuOr")[c(3,4)])+
  theme(axis.title = element_text(size=15),legend.title = element_text(size=15),axis.text = element_text(size=12),
        legend.text=element_text(size=12))
fig1

# load state and province boundaries
us1 <- gadm(country="USA",level=1)
us2 <- gadm(country="USA",level=2)
us1.sf <- st_as_sf(us1)
us2.sf <- st_as_sf(us2)
# get coordinates into correct format
us1.sf <- st_transform(us1.sf,crs=4326)
us2.sf <- st_transform(us2.sf,crs=4326)

ggplot()+
  geom_sf(data = us2.sf, color = "black", fill = NA, size = 0.3)+
  coord_sf(xlim = c(-112.4, -111.8),ylim = c(44, 44.6))+
  geom_point(data=SS.clim.data,aes(x=lon,y=lat,shape=as.factor(elevation)),size=2)+
  theme_minimal()
ggplot()+
  geom_sf(data = us1.sf, color = "black", fill = NA, size = 0.5)+
  geom_sf(data = us2.sf, color = "black", fill = NA, size = 0.1)+
  coord_sf(xlim = c(-118, -110),ylim = c(41, 49))+
  theme_minimal()+
  geom_point(aes(x=-112.222595, y=44.251839),size=6)

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
fig2

#### Fig 3: PCA showing trait dispersion -----------------

fig3 <- ggplot() +
  geom_point(data=pca_data,aes(x=PC1,y=PC2,color=elevation),alpha=0.7,size=3)+
  theme_minimal()+
  geom_segment(data=loads,aes(x=0,y=0,xend=PC1*4,yend=PC2*4),color="black",
               arrow = arrow(length=unit(0.2,"cm")))+
  geom_text(data=loads,aes(x=PC1*4.4,y=PC2*4.5,label=c("SLA","SRL","height")),color="black",size=5)+
  scale_shape_manual(values=c(2,19))+
  scale_color_manual(values=brewer.pal(11,"PuOr")[c(3,10)])+
  labs(color="Elevation")+
  theme(legend.title = element_text(size=15),legend.text = element_text(size=12),axis.title=element_text(size=15))
fig3

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
fig4

# save output
if(!dir.exists(paste0("./../outputs/", Sys.Date(),"/"))) dir.create(paste0("./../outputs/", Sys.Date(),"/"))
ggsave(fig1,file=paste0("./../outputs/", Sys.Date(),"/fig1.jpeg"),height=6,width=7)
ggsave(fig2,file = paste0("./../outputs/", Sys.Date(),"/fig2.jpeg"),height = 6,width = 10)
ggsave(fig3,file = paste0("./../outputs/", Sys.Date(),"/fig3.jpeg"),height = 5,width = 7)
ggsave(fig4, file=paste0("./../outputs/", Sys.Date(),"/fig4.jpeg"),height=5,width=7)

##### Supplemental figures --------------

### Phenology figures

# formatting data
etraits <- gdata[-which(is.na(gdata$emergence)),]
etraits$emergence <- as.numeric(strftime(etraits$emergence, format = "%V"))

ftraits <- gdata[-which(is.na(gdata$flower)),]
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
low.sm <- mean(gdata$mean.season.VWC[gdata$elevation=="low"])
high.sm <- mean(gdata$mean.season.VWC[gdata$elevation=="high"])
open.sm <- mean(gdata$mean.season.VWC[gdata$patch=="open"])
shrub.sm <- mean(gdata$mean.season.VWC[gdata$patch=="shrub"])

# calculate temperature differences
low.st <- mean(gdata$mean.season.temp[gdata$elevation=="low" & !is.na(gdata$mean.season.temp)])
high.st <- mean(gdata$mean.season.temp[gdata$elevation=="high" & !is.na(gdata$mean.season.temp)])
open.st <- mean(gdata$mean.season.temp[gdata$patch=="open" & !is.na(gdata$mean.season.temp)])
shrub.st <- mean(gdata$mean.season.temp[gdata$patch=="shrub" & !is.na(gdata$mean.season.temp)])

# two way anova on soil moisture differences between patches and elevations
mean.sm.aov <- aov(mean.season.VWC ~ elevation + patch, data=gdata)
summary(mean.sm.aov)
# two way anova on soil temperature differences between patches and elevations
mean.st.aov <- aov(mean.season.temp ~ elevation + patch, data=gdata)
summary(mean.st.aov)

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
moisture.data <- pivot_longer(data=gdata,cols = c("mean.season.VWC","min.season.VWC","max.season.VWC"),names_to = "anomoly",values_to = "percent_VWC")
temp.data <- pivot_longer(data=gdata,cols = c("mean.season.temp","min.season.temp","max.season.temp"),names_to = "anomoly",values_to = "degrees_C")
temp.data$degrees_C <- (temp.data$degrees_C - 32) * (5/9)

# plot soil moisture
moisture.plot <- ggplot(moisture.data,aes(x=as.factor(elevation),y=percent_VWC,fill=as.factor(patch)))+
  facet_grid(~anomoly,labeller=labeller(anomoly=c("max.season.VWC"="maximum","mean.season.VWC"="mean",
                                                   "min.season.VWC"="minimum")))+
  geom_boxplot()+
  labs(x="Elevation",y="Average season soil moisture (% VWC)",fill="Patch")+
  theme_minimal()+
  theme(strip.text = element_text(size=12),axis.title = element_text(size=15),axis.text = element_text(size=12),
        legend.title = element_text(size=15),legend.text=element_text(size=12))
moisture.plot

# plot soil temperature
temp.plot <- ggplot(temp.data,aes(x=as.factor(elevation),y=degrees_C,fill=as.factor(patch)))+
  facet_grid(~anomoly,labeller=labeller(anomoly=c("max.season.temp"="maximum","mean.season.temp"="mean",
                                                  "min.season.temp"="minimum")))+
  geom_boxplot()+
  labs(x="Elevation",y=expression("Average season soil temperature (" * degree * "C)"),fill="Patch")+
  theme_minimal()+
  theme(strip.text = element_text(size=12),axis.title = element_text(size=15),axis.text = element_text(size=12),
        legend.title = element_text(size=15),legend.text=element_text(size=12))
temp.plot

### Show that GDD and snowmelt timing correlated to MAT

# MAT
ss.clim <- read.csv("./../clean_data/fig1_SS_climate.csv")
ss.clim$site <- c("low_gate","low_north","low_east","high_gate","high_north","high_east")

# snowmelt
snowmelt.dates <- snow.absent %>% group_by(year,site,patch) %>% summarize(snowmelt = min(date))
snowmelt.dates$snowmelt <- yday(snowmelt.dates$snowmelt)
snowmelt.compare <- snowmelt.dates[-c(4,11),]
snowmelt.compare <- snowmelt.compare %>% group_by(site) %>% summarize(avg.snowmelt.date = mean(as.numeric(snowmelt)))

# GDD (from temp data)
t_base <- 32
gdd_df <- soil.temp %>%
  group_by(date,site) %>%
  # Average 6-hourly temp for the day, subtract base, divide by 4 intervals
  summarize(daily_gdd = max(0, mean(temp) - t_base))

gdd_df$year <- year(ymd(gdd_df$date))
gdd_df$month <- month(ymd(gdd_df$date))
gdd_df <- gdd_df[which(gdd_df$month<8),]

gdd.sum <- gdd_df %>% group_by(site,year) %>% summarize(cum.gdd = sum(daily_gdd))
gdd.sum <- gdd.sum[-2,]
gdd.sum <- gdd.sum %>% group_by(site) %>% summarize(mean.cum.gdd = mean(cum.gdd))

# putting it all together
snow.temp.gdd <- left_join(snowmelt.compare,ss.clim,by=c("site"))
snow.temp.gdd <- left_join(snow.temp.gdd,gdd.sum,by=c("site"))

cor(snow.temp.gdd$avg.snowmelt.date,snow.temp.gdd$MAT.F)
cor(snow.temp.gdd$MAT.F,snow.temp.gdd$mean.cum.gdd)

p1 <- ggplot(snow.temp.comp,aes(x=MAT.F,y=avg.snowmelt.date))+
  geom_point()+
  geom_smooth(method="lm")+
  geom_text(x=43.5,y=100,label = "r = -0.96",size=5)+
  labs(x="MAT (F)",y="Avg. snowmelt day of year")
p2 <- ggplot(snow.temp.gdd,aes(x=MAT.F,y=mean.cum.gdd))+
  geom_point()+
  geom_smooth(method="lm")+
  geom_text(x=43.5,y=1500,label = "r = 0.84",size=5)+
  labs(x="MAT (F)",y="Avg. cumulative GDD")
ggarrange(p1,p2, ncol=1,nrow=2)

# testing for intraspecific facilitation

# biomass x proportion cheatgrass
# bdata$prop.cheat <- 0
# for(i in 1:nrow(bdata)){
#   if((bdata$n.con.nb[i]+bdata$n.het.nb[i])==0){bdata$prop.cheat[i]<-0}else{
#     bdata$prop.cheat[i] <- (bdata$n.con.nb[i] / (bdata$n.con.nb[i]+bdata$n.het.nb[i]))
#   }
# }
# # plot
# ggplot(bdata[which(bdata$subplot=="C"),],aes(x=prop.cheat,y=log(aboveground_mass)))+
#   geom_point()+
#   geom_smooth(method="lm")
# # quick lm - not significant
# bdata$prop.cheat <- scale(bdata$prop.cheat)
# mod <- lmer(log(aboveground_mass) ~ prop.cheat + elevation + patch + (1|siteyear),data=bdata)
# summary(mod)
# 
# # biomass x cheatgrass abundance
# # plot
# ggplot(bdata[which(bdata$subplot=="C"),],aes(x=n.con.nb,y=log(aboveground_mass)))+
#   geom_point()+
#   geom_smooth(method="lm")
# # quick lm - not significant
# bdata$n.con.nb <- scale(bdata$n.con.nb)
# mod <- lmer(log(aboveground_mass) ~ n.con.nb + elevation + patch + (1|siteyear),data=bdata)
# summary(mod)
