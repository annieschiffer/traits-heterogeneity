## Make figures

# import data
data <- read.csv("./../clean_data/all_data_combined.csv")
data <- data[data$species=="BRTE",]

ltraits <- data[-which(is.na(data$total_leaf_area)),]
ltraits$total_leaf_area <- ltraits$total_leaf_area/100
ltraits$SLA <- as.numeric(ltraits$total_leaf_area)/as.numeric(ltraits$total_leaf_mass)

rtraits <- data[-which(is.na(data$length_cm)),]
rtraits <- rtraits[-which(is.na(rtraits$root_mass)),]
rtraits <- rtraits[-which(rtraits$root_mass==0),]
rtraits$SRL <- as.numeric(rtraits$length_cm)/as.numeric(rtraits$root_mass)

htraits <- data[-which(is.na(data$max.height)),]
htraits <- htraits[-which(htraits$max.height==0),]

etraits <- data[-which(is.na(data$emergence)),]
etraits$emergence <- as.numeric(strftime(etraits$emergence, format = "%V"))

ftraits <- data[-which(is.na(data$flower)),]
ftraits$flower <- as.numeric(strftime(ftraits$flower, format = "%V"))

##### Fig 2 --------------

# fixed effects: germination rates & biomass
get.coeff <- function(mod){
  modfe <- fixef(mod)
  modfe <- modfe[-1]
  modfe <- unname(modfe)
  se_fixed <- sqrt(diag(vcov(mod)))
  se_fixed <- se_fixed[-1]
  modcoeff <- as.data.frame(cbind(modfe,se_fixed))
  modcoeff$coeff <- rownames(modcoeff)
  rownames(modcoeff) <- NULL
  return(modcoeff)
}
germcoeff <- get.coeff(germ)
masscoeff <- get.coeff(mass)
germcoeff$vital.rate <- "Germination"
masscoeff$vital.rate <- "Biomass"

coeffs <- rbind(germcoeff,masscoeff)

rect_data1 <- data.frame(xmin = -Inf, xmax = Inf,ymin = 3.75, ymax = 5.5)
rect_data2 <- data.frame(xmin = -Inf, xmax = Inf,ymin = 2.25, ymax = 3.75)
rect_data3 <- data.frame(xmin = -Inf, xmax = Inf,ymin = 0.5, ymax = 2.25)
rect_data4 <- data.frame(xmin = -Inf, xmax = Inf,ymin = 6.5, ymax = 9)
rect_data5 <- data.frame(xmin = -Inf, xmax = Inf,ymin = 4.3, ymax = 6.5)
rect_data6 <- data.frame(xmin = -Inf, xmax = Inf,ymin = 0.5, ymax = 4.3)

germ.plot <-ggplot(coeffs[which(coeffs$vital.rate=="Germination"),],aes(x=modfe,y=coeff))+
  geom_pointrange(aes(xmin=modfe-se_fixed,xmax=modfe+se_fixed))+
  theme_minimal()+
  geom_vline(xintercept=0,lty="dashed",alpha=0.5)+
  scale_y_discrete(limits=c("nb.number:patchshrub","elevationlow:nb.number",
                            "nb.number","patchshrub","elevationlow"),
                   labels=c("Neighbors x Patch (shrub)","Neighbors x Elevation (low)",
                            "Number of neighbors","Patch (shrub)","Elevation (low)"))+
  theme(axis.title.y=element_blank(),axis.text.y = element_text(size=12,color="black"))+
  labs(x="Model estimate",title="Germination")+
  geom_rect(data = rect_data1, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = brewer.pal(11,"PuOr")[4], alpha = 0.2, inherit.aes = FALSE)+
  geom_rect(data = rect_data2, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = brewer.pal(11,"PuOr")[9], alpha = 0.2, inherit.aes = FALSE)+
  geom_rect(data = rect_data3, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = brewer.pal(11,"PuOr")[2], alpha = 0.2, inherit.aes = FALSE)+
  annotate("text",x=3.5,y=4.7,label="Environment",angle=270)+
  annotate("text",x=3.5,y=3,label="Competition",angle=270)+
  annotate("text",x=3.5,y=1.4,label="Interaction",angle=270)

biomass.plot <- ggplot(coeffs[which(coeffs$vital.rate=="Biomass"),],aes(x=modfe,y=coeff))+
  geom_pointrange(aes(xmin=modfe-se_fixed,xmax=modfe+se_fixed))+
  theme_minimal()+
  geom_vline(xintercept=0,lty="dashed",alpha=0.5)+
  scale_y_discrete(limits=c("patchshrub:n.het.nb","n.het.nb:elevationlow","n.con.nb:patchshrub",
  "n.con.nb:elevationlow","n.het.nb","n.con.nb","patchshrub","elevationlow"),
                   labels=c("Heterospecific x Patch (shrub)","Heterospecific x Elevation (low)",
                            "Conspecific x Patch (shrub)","Conspecific x Elevation (low)",
                            "Heterospecific neighbors","Conspecific neighbors",
                            "Patch (shrub)","Elevation (low)"))+
  theme(axis.title.y=element_blank(),axis.text.y = element_text(size=12,color="black"))+
  labs(x="Model estimate",title="Biomass")+
  scale_x_continuous(breaks=seq(-8,8,by=2))+
  geom_rect(data = rect_data4, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = brewer.pal(11,"PuOr")[4], alpha = 0.2, inherit.aes = FALSE)+
  geom_rect(data = rect_data5, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = brewer.pal(11,"PuOr")[9], alpha = 0.2, inherit.aes = FALSE)+
  geom_rect(data = rect_data6, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = brewer.pal(11,"PuOr")[2], alpha = 0.2, inherit.aes = FALSE)+
  annotate("text",x=8.75,y=7.8,label="Environment",angle=270)+
  annotate("text",x=8.75,y=5.35,label="Competition",angle=270)+
  annotate("text",x=8.75,y=2.4,label="Interaction",angle=270)

vital.rates <- ggarrange(germ.plot,biomass.plot,widths = c(0.7,1))

#### Fig 3 -----------------

# effect of neighbors in low (favorable) vs. high (stressful) elevation

# formatting
biomass <- data[!is.na(data$aboveground_mass),]
biomass$siteyear <- paste0(biomass$site,biomass$year)
biomass$siteyear <- as.numeric(as.factor(biomass$siteyear))

biomass <- pivot_longer(biomass,cols = c("n.het.nb","n.con.nb"),names_to = "Neighbors",
                        values_to = "Neighbor_number")

# all neighbors together
nb.plot <- ggplot(biomass,aes(x=Neighbor_number,y=log(aboveground_mass),color=elevation,shape=Neighbors))+
  geom_point(alpha=0.4)+
  geom_smooth(method="lm",aes(lty=Neighbors),se=FALSE)+
  labs(x="Number of neighbors",y="log(biomass)",color="Elevation",title=
         "Effect of neighbors on aboveground biomass")+
  scale_shape_discrete(limits=c("n.con.nb","n.het.nb"),labels=c("conspecific","heterospecific"))+
  scale_linetype_discrete(limits=c("n.con.nb","n.het.nb"),labels=c("conspecific","heterospecific"))+
  theme_classic()+
  scale_color_manual(limits=c("high","low"),values=brewer.pal(11,"PuOr")[c(9,3)])+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=12),legend.title = element_text(size=15),
        legend.text = element_text(size=12))

##### Fig 4 ----------------

# trait shifts in low vs. high and open vs. shrub

lsum <- ltraits[,c("patch","elevation","SLA")]
lsum <- rename(lsum, "value"="SLA")
lsum$trait <- "SLA"
rsum <- rtraits[,c("patch","elevation","SRL")]
rsum <- rename(rsum, "value"="SRL")
rsum$trait <- "SRL"
hsum <- htraits[,c("patch","elevation","max.height")]
hsum <- rename(hsum, "value"="max.height")
hsum$trait <- "height"

trait.sum <- rbind(lsum,rsum,hsum)
trait.sum$elevation <- factor(trait.sum$elevation, levels = c("low","high"))

trait.patch <- ggplot(trait.sum, aes(x=patch,y=log(value),fill=trait))+
  geom_violin(position = position_dodge(width=0.2))+
  theme_minimal()+
  scale_fill_manual(limits=c("SRL","SLA","height"),values=brewer.pal(11,"PuOr")[c(2,8,10)])+
  annotate("text",label="*",x=1.5,y=1,size=12)+
  labs(title = "Trait shifts across patches and elevation",x="Patch",y="log (Value)",fill="Trait")+
  theme(axis.text=element_text(size=12),axis.title = element_text(size=15),legend.position="none",
        title = element_text(size=15))
trait.elev <- ggplot(trait.sum, aes(x=elevation,y=log(value),fill=trait))+
  geom_violin(position = position_dodge(width=0.2))+
  theme_minimal()+
  scale_fill_manual(limits=c("SRL","SLA","height"),values=brewer.pal(11,"PuOr")[c(2,8,10)])+
  annotate("text",label="*",x=1.5,y=1,size=12)+
  annotate("text",label="*",x=1.5,y=5,size=12)+
  annotate("text",label="*",x=1.5,y=9.5,size=12)+
  labs(title = " ",x="Elevation",y="log (Value)",fill="Trait")+
  theme(axis.text.x=element_text(size=12),axis.title.x = element_text(size=15),axis.text.y=element_blank(),
        axis.title.y=element_blank(),legend.title = element_text(size=15),
        legend.text = element_text(size=12),title=element_text(size=15))

trait.shifts <- ggarrange(trait.patch,trait.elev,widths=c(0.8,1))

#### Fig 5 -----------------

# scatterplot for number of neighbors x continuous trait variables
format.traits <- function(trait.data,col,name){
  d <- trait.data[,c(name,"n.het.nb","n.con.nb")]
  d$trait <- name
  d$value <- col
  d$value <- scale(d$value)
  d <- d[,-1]
  return(d)
}

lnb <- format.traits(ltraits,ltraits$SLA,"SLA")
rnb <- format.traits(rtraits,rtraits$SRL,"SRL")
hnb <- format.traits(htraits,htraits$max.height,"max.height")

trait.nb <- rbind(lnb,rnb,hnb)
trait.nb <- pivot_longer(trait.nb,cols=c("n.het.nb","n.con.nb"),names_to = "Neighbors",values_to = "n.nb")

trait.nb <- ggplot(trait.nb,aes(x=n.nb,y=value,color=trait,shape=Neighbors))+
  geom_point(alpha=0.4)+
  geom_smooth(method="lm",aes(lty=Neighbors),se=FALSE)+
  scale_color_manual(limits=c("SLA","SRL","max.height"),values=c(brewer.pal(11,"PuOr")[c(2,8,10)]),
                     labels=c("SLA","SRL","height"))+
  labs(x="Number of neighbors",y="Scaled trait value",color="Trait",title="Effect of neighbors on traits")+
  scale_shape_discrete(limits=c("n.con.nb","n.het.nb"),labels=c("conspecific","heterospecific"))+
  scale_linetype_discrete(limits=c("n.con.nb","n.het.nb"),labels=c("conspecific","heterospecific"))+
  theme_classic()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=12),legend.title = element_text(size=15),
        legend.text = element_text(size=12))

# save output
if(!dir.exists(paste0("./../results/", Sys.Date(),"/"))) dir.create(paste0("./../results/", Sys.Date(),"/"))
ggsave(vital.rates,file = paste0("./../results/", Sys.Date(),"/vital_estimates.jpeg"),height = 4,width = 10)
ggsave(nb.plot,file = paste0("./../results/",Sys.Date(),"/SGH_plot.jpeg"),height = 5,width = 6)
ggsave(trait.shifts, file = paste0("./../results/",Sys.Date(),"/trait_shifts.jpeg"),height=5,width=8)
ggsave(trait.nb,file=paste0("./../results/",Sys.Date(),"/traits_neighbors.jpeg"),height=5,width=6)


##### Supplemental figures --------------

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
