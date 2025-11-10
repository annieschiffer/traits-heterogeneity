### fitting models with traits as response

data <- read.csv("./../clean_data/all_data_combined.csv")
data <- data[data$species=="BRTE",]

## SLA

# calculate SLA
ltraits <- data[-which(is.na(data$total_leaf_area)),]
ltraits$total_leaf_area <- ltraits$total_leaf_area/100
ltraits$SLA <- as.numeric(ltraits$total_leaf_area)/as.numeric(ltraits$total_leaf_mass)
# set up siteyear random effect & scaled predictors
ltraits$siteyear <- paste0(ltraits$site,ltraits$year)
ltraits$siteyear <- as.numeric(as.factor(ltraits$siteyear))
ltraits$n.con.nb <- as.numeric(scale(ltraits$n.con.nb))
ltraits$n.het.nb <- as.numeric(scale(ltraits$n.het.nb))

# look at distribution of trait
hist(log(ltraits$SLA)) # log-normal
# model
lfit <- lmer(log(SLA) ~ n.het.nb + n.con.nb + patch + elevation + (1|siteyear),data=ltraits)
summary(lfit)
res <- simulateResiduals(lfit)
base::plot(res) # residuals good enough

## SRL

# calculate SRL
rtraits <- data[-which(is.na(data$length_cm)),]
rtraits <- rtraits[-which(is.na(rtraits$root_mass)),]
rtraits <- rtraits[-which(rtraits$root_mass==0),]
rtraits$SRL <- as.numeric(rtraits$length_cm)/as.numeric(rtraits$root_mass)

# set up siteyear random effect & scaled predictors
rtraits$siteyear <- paste0(rtraits$site,rtraits$year)
rtraits$siteyear <- as.numeric(as.factor(rtraits$siteyear))
rtraits$n.con.nb <- as.numeric(scale(rtraits$n.con.nb))
rtraits$n.het.nb <- as.numeric(scale(rtraits$n.het.nb))

# look at distribution of trait
hist(log(rtraits$SRL)) # log-normal
# fit model
rfit <- lmer(log(SRL) ~ n.het.nb + n.con.nb + patch + elevation + (1|siteyear),data=rtraits)
summary(rfit)
res <- simulateResiduals(rfit)
base::plot(res) # residuals good


## height

htraits <- data[-which(is.na(data$max.height)),]
htraits <- htraits[-which(htraits$max.height==0),]

# set up siteyear random effect & scaled predictors
htraits$siteyear <- paste0(htraits$site,htraits$year)
htraits$siteyear <- as.numeric(as.factor(htraits$siteyear))
htraits$n.con.nb <- as.numeric(scale(htraits$n.con.nb))
htraits$n.het.nb <- as.numeric(scale(htraits$n.het.nb))

# look at distribution of trait
hist(log(htraits$max.height)) # log-normal
# fit model
hfit <-lmer(log(max.height) ~ n.het.nb + n.con.nb + patch + elevation + (1|siteyear),data=htraits)
summary(hfit)
res <- simulateResiduals(hfit)
base::plot(res) # residuals good


## emergence phenology

etraits <- data[-which(is.na(data$emergence)),]
etraits$emergence <- as.numeric(strftime(etraits$emergence, format = "%V"))

# set up siteyear random effect & scaled predictors
etraits$siteyear <- paste0(etraits$site,etraits$year)
etraits$siteyear <- as.numeric(as.factor(etraits$siteyear))
etraits$n.con.nb <- as.numeric(scale(etraits$n.con.nb))
etraits$n.het.nb <- as.numeric(scale(etraits$n.het.nb))

# look at distribution of trait
hist(etraits$emergence) # this is causing issues - come back to this


## flowering phenology

ftraits <- data[-which(is.na(data$flower)),]
ftraits$flower <- as.numeric(strftime(ftraits$flower, format = "%V"))

# set up siteyear random effect & scaled predictors
ftraits$siteyear <- paste0(ftraits$site,ftraits$year)
ftraits$siteyear <- as.numeric(as.factor(ftraits$siteyear))
ftraits$n.con.nb <- as.numeric(scale(ftraits$n.con.nb))
ftraits$n.het.nb <- as.numeric(scale(ftraits$n.het.nb))

# look at distribution of trait
hist(ftraits$flower) # maybe normal? check residuals

ffit <-lmer(flower ~ n.het.nb + n.con.nb + patch + elevation + (1|siteyear),data=ftraits)
res <- simulateResiduals(ffit)
base::plot(res) # residuals bad

