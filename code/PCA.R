
# import data
data <- read.csv("./../clean_data/all_data_combined.csv")
data <- data[data$species=="BRTE",]

# only complete data set
bdata <- data[!is.na(data$aboveground_mass),]
traits <- bdata[,c("max.height","length_cm","root_mass","total_leaf_area","total_leaf_mass")]
bdata.complete <- bdata[which(complete.cases(traits)==TRUE),]
#bdata.complete <- bdata.complete[bdata.complete$subplot=="R",]

# calculate SLA
bdata.complete$total_leaf_area <- bdata.complete$total_leaf_area/100
bdata.complete$SLA <- as.numeric(bdata.complete$total_leaf_area)/as.numeric(bdata.complete$total_leaf_mass)

# calculate SRL
bdata.complete$SRL <- as.numeric(bdata.complete$length_cm)/as.numeric(bdata.complete$root_mass)
bdata.complete <- bdata.complete[-which(bdata.complete$SRL==Inf),]

# just get traits
bdata.complete$max.height <- as.numeric(bdata.complete$max.height)
pca.data <- bdata.complete[,c("SLA","SRL","max.height")]

pca <- prcomp(pca.data,center=TRUE,scale. = TRUE)

summary(pca)
plot(pca)

# add pc scores to data frame
bdata.complete$PC1 <- pca$x[,1]
bdata.complete$PC2 <- pca$x[,2]

write.csv(bdata.complete,"./../clean_data/pca_data.csv",row.names = FALSE)

# extract loadings
loads <- as.data.frame(pca$rotation[,1:2])

# testing for differences in dispersion
elev.disp <- betadisper(dist(bdata.complete[,c("PC1","PC2")]), group = bdata.complete$elevation)
anova(elev.disp)
# functional dispersion at low > at high

# pat.disp <- betadisper(dist(bdata.complete[,c("PC1","PC2")]), group = bdata.complete$patch)
# anova(pat.disp) # nothing significant

rm(bdata.complete)
