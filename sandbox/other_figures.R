## Abandoned figures and analyses

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

# testing for intraspecific facilitation

# biomass x proportion cheatgrass
bdata$prop.cheat <- 0
for(i in 1:nrow(bdata)){
  if((bdata$n.con.nb[i]+bdata$n.het.nb[i])==0){bdata$prop.cheat[i]<-0}else{
    bdata$prop.cheat[i] <- (bdata$n.con.nb[i] / (bdata$n.con.nb[i]+bdata$n.het.nb[i]))
  }
}
# plot
ggplot(bdata[which(bdata$subplot=="C"),],aes(x=prop.cheat,y=log(aboveground_mass)))+
  geom_point()+
  geom_smooth(method="lm")
# quick lm - not significant
bdata$prop.cheat <- scale(bdata$prop.cheat)
mod <- lmer(log(aboveground_mass) ~ prop.cheat + elevation + patch + (1|siteyear),data=bdata)
summary(mod)

# biomass x cheatgrass abundance
# plot
ggplot(bdata[which(bdata$subplot=="C"),],aes(x=n.con.nb,y=log(aboveground_mass)))+
  geom_point()+
  geom_smooth(method="lm")
# quick lm - not significant
bdata$n.con.nb <- scale(bdata$n.con.nb)
mod <- lmer(log(aboveground_mass) ~ n.con.nb + elevation + patch + (1|siteyear),data=bdata)
summary(mod)
