## Make figures

##### Fig 2 --------------
load("./../outputs/2026-02-25/germ.fit.rda")
load("./../outputs/2026-02-25/mass.fit.rda")

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
  labs(x="Scaled effect size",title="Germination")+
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

# import data
all.data <- read.csv("./../clean_data/all_data_combined.csv")
all.data <- all.data[all.data$species=="BRTE",]

germ.rates <- all.data %>% group_by(elevation,patch,subplot) %>% summarize(germ.rate = sum(germination)/n())
sgh.elev <- ggplot(germ.rates,aes(x=elevation,y=germ.rate,fill=subplot))+
  geom_bar(stat = "identity",position=position_dodge())+
  labs(x="Environment",y="Germination rate",fill="Neighbors",
       title="Germination rates across stress gradient")+
  scale_fill_manual(limits=c("R","C"),
                    labels=c("absent","present"),
                    values=brewer.pal(11,"PuOr")[c(2,10)])+
  scale_x_discrete(limits=c("low","high"),
                   labels=c("favorable","stressful"))+
  theme_minimal()+
  theme(axis.title = element_text(size=15),axis.text=element_text(size=12),title=element_text(size=15,face="bold"),
        legend.title=element_text(size=15),legend.text=element_text(size=12))
sgh.elev

##### Fig 4 ----------------

# trait shifts in low vs. high and open vs. shrub - only SRL varied with anything

rtraits <- all.data[-which(is.na(all.data$length_cm)),]
rtraits <- rtraits[-which(is.na(rtraits$root_mass)),]
rtraits <- rtraits[-which(rtraits$root_mass==0),]
rtraits$SRL <- as.numeric(rtraits$length_cm)/as.numeric(rtraits$root_mass)

rsum <- rtraits[,c("patch","elevation","n.het.nb","n.con.nb","SRL")]
rsum<- pivot_longer(rsum,cols = c(n.het.nb,n.con.nb),names_to = "Neighbors",values_to = "Number")
rsum <- rename(rsum, "value"="SRL")
rsum$trait <- "SRL"
rsum$elevation <- factor(rsum$elevation, levels = c("low","high"))

trait.scatter <- ggplot(rsum,aes(x=Number,y=log(value),color=Neighbors,linetype=patch))+
  geom_point()+
  geom_smooth(method="lm",se=FALSE)+
  scale_color_manual(limits=c("n.con.nb","n.het.nb"),
                     labels=c("conspecific","heterospecific"),
    values=brewer.pal(11,"PuOr")[c(3,9)])+
  scale_linetype_manual(limits=c("open","shrub"),values=c("dashed","solid"))+
  theme_minimal()+
  labs(x="Number of neighbors",y="log(SRL)",title="Effect of patch and neighbors on SRL")+
  theme(legend.title = element_blank(),legend.text = element_text(size=12),axis.text=element_text(size=12),axis.title=element_text(size=15),
        title=element_text(size=15,face="bold"))
trait.scatter

#### Fig 5? -----------------




# save output
if(!dir.exists(paste0("./../outputs/", Sys.Date(),"/"))) dir.create(paste0("./../outputs/", Sys.Date(),"/"))
ggsave(vital.rates,file = paste0("./../outputs/", Sys.Date(),"/vital_estimates.jpeg"),height = 6,width = 10)
ggsave(sgh.elev,file = paste0("./../outputs/",Sys.Date(),"/SGH_plot.jpeg"),height = 5,width = 6)
ggsave(trait.scatter, file = paste0("./../outputs/",Sys.Date(),"/SRL_scatter.jpeg"),height=5,width=6)



##### Supplemental figures --------------

etraits <- data[-which(is.na(data$emergence)),]
etraits$emergence <- as.numeric(strftime(etraits$emergence, format = "%V"))

ftraits <- data[-which(is.na(data$flower)),]
ftraits$flower <- as.numeric(strftime(ftraits$flower, format = "%V"))

ecount <- etraits %>% group_by(patch,elevation,subplot,emergence) %>% summarize(emerged = n())
fcount <- ftraits %>% group_by(patch,elevation,subplot,flower) %>% summarise(flowered=n())

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

ggarrange(eplot,fplot,nrow=2,ncol=1)

# Elevation & environment

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

hist(snow.absent$date[snow.absent$site=="low_north" & snow.absent$year==2025],breaks="days")

snowmelt.dates <- snow.absent %>% group_by(year,elevation,patch) %>% summarize(snowmelt = min(date))
snowmelt.dates$snowmelt <- yday(snowmelt.dates$snowmelt)

low.date <- mean(snowmelt.dates$snowmelt[snowmelt.dates$elevation=="low"])
high.date <- mean(snowmelt.dates$snowmelt[snowmelt.dates$elevation=="high"])
open.date <- mean(snowmelt.dates$snowmelt[snowmelt.dates$patch=="open"])
shrub.date <- mean(snowmelt.dates$snowmelt[snowmelt.dates$patch=="shrub"])

# ggplot(data,aes(x=as.factor(year),y=mean.season.VWC,fill=as.factor(elevation)))+
#   geom_boxplot()
# ggplot(data,aes(x=as.factor(year),y=min.season.VWC,fill=as.factor(elevation)))+
#   geom_boxplot()
# ggplot(data,aes(x=as.factor(year),y=max.season.VWC,fill=as.factor(elevation)))+
#   geom_boxplot()
# 
# ggplot(data,aes(x=as.factor(elevation),y=mean.season.temp))+
#   geom_boxplot()
# ggplot(data,aes(x=as.factor(elevation),y=min.season.temp))+
#   geom_boxplot()
# ggplot(data,aes(x=as.factor(elevation),y=max.season.temp))+
#   geom_boxplot()
# 
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
