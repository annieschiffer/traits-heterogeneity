## Make figures

# # import data
# data <- read.csv("./../clean_data/all_data_combined.csv")
# data <- data[data$species=="BRTE",]
# data$nb.number <- data$n.con.nb + data$n.het.nb

## Fig 1
# multi-panel pictures of experimental setup

## Fig 2 
# fixed effects: germination rates & biomass (see drawing)
get.coeff <- function(mod){
  modfe <- fixef(mod)
  modfe <- modfe[-1]
  modfe <- unname(modfe)
  modci <- confint(mod, level=0.95)
  if(class(mod)=="glmerMod"){modconf<-as.data.frame(modci[-c(1,2),])}else{modconf<-as.data.frame(modci[-c(1:3),])}
  modconf$coeff <- row.names(modconf)
  rownames(modconf) <- NULL
  modcoeff <- cbind(modfe,modconf)
  return(modcoeff)
}
germcoeff <- get.coeff(germ)
masscoeff <- get.coeff(mass)

germcoeff$vital.rate <- "Germination"
masscoeff$vital.rate <- "Biomass"
coeffs <- rbind(germcoeff,masscoeff)

rect_data1 <- data.frame(
  xmin = -Inf, xmax = Inf,
  ymin = 3.75, ymax = 5.5 # Extend across entire y-axis
)
rect_data2 <- data.frame(
  xmin = -Inf, xmax = Inf,
  ymin = 2.25, ymax = 3.75 # Extend across entire y-axis
)
rect_data3 <- data.frame(
  xmin = -Inf, xmax = Inf,
  ymin = 0.5, ymax = 2.25 # Extend across entire y-axis
)
ggplot(coeffs[which(coeffs$vital.rate=="Germination"),],aes(x=modfe,y=coeff))+
  geom_pointrange(aes(xmin=`2.5 %`,xmax=`97.5 %`))+
  theme_minimal()+
  geom_vline(xintercept=0,lty="dashed",alpha=0.5)+
  scale_y_discrete(limits=c("nb.number:patchshrub","elevationlow:nb.number",
                            "nb.number","patchshrub","elevationlow"),
                   labels=c("Neighbors x Patch","Neighbors x Elevation",
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

rect_data4 <- data.frame(
  xmin = -Inf, xmax = Inf,
  ymin = 6.5, ymax = 9 # Extend across entire y-axis
)
rect_data5 <- data.frame(
  xmin = -Inf, xmax = Inf,
  ymin = 4.3, ymax = 6.5 # Extend across entire y-axis
)
rect_data6 <- data.frame(
  xmin = -Inf, xmax = Inf,
  ymin = 0.5, ymax = 4.3 # Extend across entire y-axis
)
ggplot(coeffs[which(coeffs$vital.rate=="Biomass"),],aes(x=modfe,y=coeff))+
  geom_pointrange(aes(xmin=`2.5 %`,xmax=`97.5 %`))+
  theme_minimal()+
  geom_vline(xintercept=0,lty="dashed",alpha=0.5)+
  scale_y_discrete(limits=c("patchshrub:n.het.nb","n.het.nb:elevationlow","n.con.nb:patchshrub",
  "n.con.nb:elevationlow","n.het.nb","n.con.nb","patchshrub","elevationlow"),
                   labels=c("Heterospecific neighbors x Patch","Heterospecific neighbors x Elevation",
                            "Conspecific neighbors x Patch","Conspecific neighbors x Elevation",
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


## Fig 3
# biomass reveals stress gradient hypothesis
# effect of conspecifics in low (favorable) vs. high (stressful) elevation
# import data
data <- read.csv("./../clean_data/all_data_combined.csv")
data <- data[data$species=="BRTE",]

# formatting
biomass <- data[!is.na(data$aboveground_mass),]
biomass$siteyear <- paste0(biomass$site,biomass$year)
biomass$siteyear <- as.numeric(as.factor(biomass$siteyear))

biomass <- pivot_longer(biomass,cols = c("n.het.nb","n.con.nb"),names_to = "Neighbors",
                        values_to = "Neighbor_number")

# all neighbors together
ggplot(biomass,aes(x=Neighbor_number,y=log(aboveground_mass),color=elevation,shape=Neighbors))+
  geom_point(alpha=0.4)+
  geom_smooth(method="lm",aes(lty=Neighbors),se=FALSE)+
  labs(x="Number of neighbors",y="log(biomass)",color="Elevation")+
  scale_shape_discrete(limits=c("n.con.nb","n.het.nb"),labels=c("conspecific","heterospecific"))+
  scale_linetype_discrete(limits=c("n.con.nb","n.het.nb"),labels=c("conspecific","heterospecific"))+
  theme_classic()+
  scale_color_manual(limits=c("high","low"),values=brewer.pal(11,"PuOr")[c(9,3)])

# color=brewer.pal(11,"PuOr")[c(1,3,4,8,10)]

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
