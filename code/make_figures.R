## Make figures

#### load data and custom functions
source("fig_functions.R")
load("./../outputs/2026-03-04/stan_fits/germ.fit.rda")
load("./../outputs/2026-03-04/stan_fits/mass.fit.rda")
load("./../outputs/2026-03-04/stan_fits/TbyE.mass.fit.rda")
# import data
gdata <- read.csv("./../clean_data/all_data_combined.csv")
gdata <- gdata[gdata$species=="BRTE",]
all.data <- gdata
bdata.complete <- read.csv("./../clean_data/complete_data_biomass.csv")

##### Fig 2 --------------

# for binary variables, first alphabetical level = reference
# output shows difference moving from high to low elevation, difference moving from open to shrub patch

rect_data1 <- data.frame(xmin = -Inf, xmax = Inf,ymin = 6.5, ymax = 8.5)
rect_data2 <- data.frame(xmin = -Inf, xmax = Inf,ymin = 4.5, ymax = 6.5)
rect_data3 <- data.frame(xmin = -Inf, xmax = Inf,ymin = 0, ymax = 4.5)
rect_data4 <- data.frame(xmin = -Inf, xmax = Inf,ymin = 6.5, ymax = 8.5)
rect_data5 <- data.frame(xmin = -Inf, xmax = Inf,ymin = 4.5, ymax = 6.5)
rect_data6 <- data.frame(xmin = -Inf, xmax = Inf,ymin = 0, ymax = 4.5)

germ.results <- as.matrix(germ.fit,pars=c("beta"))
mass.results <- as.matrix(mass.fit,pars=c("beta"))
color_scheme_set(scheme = "darkgray")

germ.fig <- mcmc_intervals(germ.results) +
  scale_y_discrete(limit=rev(c("beta[1]","beta[2]","beta[3]","beta[4]","beta[5]","beta[6]","beta[7]","beta[8]")),
                   label=rev(c("elevation (low)","patch (shrub)","conspecific neighbors","heterospecific neighbors","elevation (low) x conspecific","elevation (low) x heterospecific",
                           "patch (shrub) x conspecific","patch (shrub) x heterospecific")))+
  labs(x="Scaled effect size",title="Emergence")+
  theme_classic()+
  geom_rect(data = rect_data1, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = brewer.pal(11,"PuOr")[4], alpha = 0.2, inherit.aes = FALSE)+
  geom_rect(data = rect_data2, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = brewer.pal(11,"PuOr")[9], alpha = 0.2, inherit.aes = FALSE)+
  geom_rect(data = rect_data3, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = brewer.pal(11,"PuOr")[2], alpha = 0.2, inherit.aes = FALSE)+
  annotate("text",x=11.35,y=7.5,label="Environment",angle=270)+
  annotate("text",x=11.35,y=5.5,label="Competition",angle=270)+
  annotate("text",x=11.35,y=2.5,label="Interaction",angle=270)+
  theme(axis.text = element_text(size=12),axis.title=element_text(size=12),title = element_text(size=15,face="bold"))
germ.fig
mass.fig <- mcmc_intervals(mass.results) +
  scale_y_discrete(limit=rev(c("beta[1]","beta[2]","beta[3]","beta[4]","beta[5]","beta[6]","beta[7]","beta[8]")),
                   label=rev(c("elevation (low)","patch (shrub)","conspecific neighbors","heterospecific neighbors","elevation (low) x conspecific","elevation (low) x heterospecific",
                               "patch (shrub) x conspecific","patch (shrub) x heterospecific")))+
  labs(x="Scaled effect size",title="Aboveground biomass")+
  theme_classic()+
  geom_rect(data = rect_data4, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = brewer.pal(11,"PuOr")[4], alpha = 0.2, inherit.aes = FALSE)+
  geom_rect(data = rect_data5, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = brewer.pal(11,"PuOr")[9], alpha = 0.2, inherit.aes = FALSE)+
  geom_rect(data = rect_data6, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = brewer.pal(11,"PuOr")[2], alpha = 0.2, inherit.aes = FALSE)+
  annotate("text",x=4.25,y=7.5,label="Environment",angle=270)+
  annotate("text",x=4.25,y=5.5,label="Competition",angle=270)+
  annotate("text",x=4.25,y=2.5,label="Interaction",angle=270)+
  theme(axis.text=element_text(size=12),axis.title = element_text(size=12),title=element_text(size=15,face="bold"),axis.text.y=element_blank())
mass.fig

vital.rates <- ggarrange(germ.fig,mass.fig,widths = c(1,0.7))
vital.rates

#### Fig 3 -----------------

# effect of neighbors in low (favorable) vs. high (stressful) elevation

# predicting emergence rates over new neighbor abundance data
predicted.data <- predict.emerg(germ.fit,gdata)

sgh <- ggplot(predicted.data,aes(x = num.nb, y = pmed,color=factor(elevation))) +
  facet_grid(~neighbors)+
  geom_line(size = 2) +
  geom_ribbon(aes(ymin = plower, ymax = pupper,fill=factor(elevation)),
              alpha = 0.1,colour=NA,show.legend = FALSE) +
  #ylim(0, 1)+
  scale_fill_manual(limits=c("0","1"),labels=c("high","low"),
                    values=brewer.pal(11,"PuOr")[c(10,3)])+
  scale_color_manual(limits=c("0","1"),labels=c("high (stressful)","low (favorable)"),
                    values=brewer.pal(11,"PuOr")[c(10,3)])+
  labs(x="Neighbor abundance",y="logit (emergence probability)",color="Elevation",title="Predicted emergence probability along stress gradient")+
  theme_minimal()+
  theme(strip.text = element_text(size=12),axis.title = element_text(size=15),axis.text=element_text(size=12),
        legend.title=element_text(size=15),legend.text=element_text(size=12),title=element_text(size=15))
sgh
##### Fig 4 ----------------

# trait shifts in low vs. high and open vs. shrub - only SRL varied with anything

rtraits <- all.data[-which(is.na(all.data$length_cm)),]
rtraits <- rtraits[-which(is.na(rtraits$root_mass)),]
rtraits <- rtraits[-which(rtraits$root_mass==0),]
rtraits$SRL <- as.numeric(rtraits$length_cm)/as.numeric(rtraits$root_mass)

rsum <- rtraits[,c("patch","elevation","n.het.nb","n.con.nb","SRL")]
srlpred <- predict.srl(SRL.fit,rsum)

rsum <- pivot_longer(rsum,cols = c("n.het.nb","n.con.nb"),names_to = "Neighbors",values_to = "num.nb")
srlpred$patch[srlpred$patch==0] <- "open"
srlpred$patch[srlpred$patch==1] <- "shrub"
srlpred <- srlpred %>% rename("Neighbors"="neighbors")

# boxplot
r.nb.patch <- rtraits[,c("patch","elevation","subplot","n.het.nb","n.con.nb","SRL")]
r.nb.patch <- r.nb.patch %>% mutate(neighborhood = case_when(n.het.nb > 0 & n.con.nb == 0 ~ "heterospecific",
                                              n.con.nb > 0 & n.het.nb == 0 ~ "conspecific",
                                              n.con.nb ==0 & n.het.nb == 0 ~ "none",
                                              n.het.nb > 0 & n.con.nb > 0 ~ "both"))
trait.box <- ggplot(r.nb.patch,aes(x=patch,y=log(SRL),fill=neighborhood)) +
  geom_boxplot()+
  scale_fill_manual(limits=c("both","conspecific","none"),
                    labels=c("both","conspecific only","none"),
                    values=brewer.pal(11,"PuOr")[c(4,8,9)])+
  labs(title="Effect of patch and neighbors on SRL")+
  theme_classic()+
  theme(axis.title.x = element_blank(),legend.title = element_blank(),axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=15),axis.text.y=element_text(size=12),legend.text = element_text(size=15),
        title=element_text(size=15,face="bold"))

#### Fig 5 -----------------

# effect of SRL on biomass across patches and elevations (scatterplot)

rtraits <- all.data[-which(is.na(all.data$length_cm)),]
rtraits <- rtraits[-which(is.na(rtraits$root_mass)),]
rtraits <- rtraits[-which(rtraits$root_mass==0),]
rtraits$SRL <- as.numeric(rtraits$length_cm)/as.numeric(rtraits$root_mass)

pred.y.srl <- predict.TbyE.biomass(TbyE.mass.fit,bdata.complete)

# scatterplot
trait.scatter <- ggplot(rtraits,aes(x=SRL,y=aboveground_mass,color=factor(elevation),shape=factor(patch)))+
  geom_line(data=pred.y.srl,aes(x=as.numeric(srl.raw),y=as.numeric(ymed),linetype=as.factor(patch)))+
  geom_point()+
  #facet_grid(~patch)+
  # scale_color_manual(limits=c("low","high"),
  #                    labels=c("low","high"),
  #                    values=brewer.pal(11,"PuOr")[c(3,9)])+
  # scale_shape_manual(limits=c("open","shrub"),values=c(1,19))+
  theme_minimal()+
  labs(x="Specific root length (logged)",y="log(biomass)",title="Effect of SRL x environment on biomass")+
  theme(legend.title = element_blank(),legend.text = element_text(size=12),axis.text=element_text(size=12),axis.title=element_text(size=15),
        title=element_text(size=15,face="bold"))
trait.scatter

# biomass & density
plot(all.data$n.con.nb,log(all.data$aboveground_mass))
plot(all.data$n.het.nb,log(all.data$aboveground_mass))

# save output
if(!dir.exists(paste0("./../outputs/", Sys.Date(),"/"))) dir.create(paste0("./../outputs/", Sys.Date(),"/"))
ggsave(vital.rates,file = paste0("./../outputs/", Sys.Date(),"/vital_estimates.jpeg"),height = 6,width = 10)
ggsave(sgh,file = paste0("./../outputs/",Sys.Date(),"/SGH_plot.jpeg"),height = 6,width = 10)
#ggsave(trait.scatter, file = paste0("./../outputs/",Sys.Date(),"/SRL_scatter.jpeg"),height=5,width=6)
ggsave(trait.box, file=paste0("./../outputs/",Sys.Date(),"/SRL_boxplot.jpeg"),height=5,width=7)

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

pheno.plot <- ggarrange(eplot,fplot,nrow=2,ncol=1)

# save plot
ggsave(pheno.plot,file = paste0("./../outputs/", Sys.Date(),"/supp_phenology.jpeg"),height = 6,width = 10)

### Soil conditions across elevation and patches

# calculate moisture differences
low.sm <- mean(data$mean.season.VWC[data$elevation=="low"]) 
high.sm <- mean(data$mean.season.VWC[data$elevation=="high"])
open.sm <- mean(data$mean.season.VWC[data$patch=="open"]) 
shrub.sm <- mean(data$mean.season.VWC[data$patch=="shrub"])

# calculate snowmelt date differences
spring.temps <- soil.temp[(soil.temp$date > "2024-03-01" & soil.temp$date < "2024-04-10") |
                         (soil.temp$date > "2025-03-01" & soil.temp$date < "2025-04-10"),]
snow.present <- spring.temps[which(spring.temps$temp < 35 & spring.temps$temp > 28),]
snow.absent <- spring.temps[-which(spring.temps$temp < 35 & spring.temps$temp > 28),]

#hist(snow.absent$date[snow.absent$site=="low_north" & snow.absent$year==2025],breaks="days")

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
  facet_grid(~anomoly)+
  geom_boxplot()+
  labs(x="Elevation",y="% VWC",fill="Patch")

temp.plot <- ggplot(temp.data,aes(x=as.factor(elevation),y=degrees_C,fill=as.factor(patch)))+
  facet_grid(~anomoly)+
  geom_boxplot()+
  labs(x="Elevation",y="degrees C",fill="Patch")

# save output
ggsave(moisture.plot,file=paste0("./../outputs/",Sys.Date(),"/supp_soil_moisture.jpeg"))
ggsave(temp.plot,file=paste0("./../outputs/",Sys.Date(),"/supp_soil_temp.jpeg"))


# # scatterplot for number of neighbors x continuous trait variables
# format.traits <- function(trait.data,col,name){
#   d <- trait.data[,c(name,"n.het.nb","n.con.nb")]
#   d$trait <- name
#   d$value <- col
#   d$value <- scale(d$value)
#   d <- d[,-1]
#   return(d)
# }
# 
# lnb <- format.traits(ltraits,ltraits$SLA,"SLA")
# rnb <- format.traits(rtraits,rtraits$SRL,"SRL")
# hnb <- format.traits(htraits,htraits$max.height,"max.height")
# 
# trait.nb <- rbind(lnb,rnb,hnb)
# htraits <- pivot_longer(htraits,cols=c("n.het.nb","n.con.nb"),names_to = "Neighbors",values_to = "n.nb")
# 
# trait.nb.plot <- ggplot(htraits,aes(x=n.nb,y=max.height,color=Neighbors,linetype=Neighbors))+
#   geom_point(alpha=0.4)+
#   geom_smooth(method="lm",se=FALSE)+
#   scale_color_manual(limits=c("n.het.nb","n.con.nb"),values=c(brewer.pal(11,"PuOr")[c(3,11)]),
#                      labels=c("heterospecific","conspecific"))+
#   labs(x="Number of neighbors",y="Maximum height (cm)",
#        title="Effect of neighbors on maximum height")+
#   #scale_shape_discrete(limits=c("n.con.nb","n.het.nb"),labels=c("conspecific","heterospecific"))+
#   scale_linetype_discrete(limits=c("n.con.nb","n.het.nb"),labels=c("conspecific","heterospecific"),
#                           guide="none")+
#   theme_classic()+
#   theme(axis.title = element_text(size=15),axis.text = element_text(size=12),legend.title = element_blank(),
#         legend.text = element_text(size=12))
# trait.nb.plot
