## Fit emergence and flowering phenology data with linear and quantile regressions

# import data
data <- read.csv("./../clean_data/all_data_combined.csv")
data <- data[data$species=="BRTE",]

## emergence phenology

etraits <- data[-which(is.na(data$emergence)),]
etraits$emergence <- as.numeric(strftime(etraits$emergence, format = "%V"))

# set up siteyear random effect & scaled predictors
etraits$siteyear <- paste0(etraits$site,etraits$year)
etraits$siteyear <- as.numeric(as.factor(etraits$siteyear))
etraits$n.con.nb <- as.numeric(scale(etraits$n.con.nb))
etraits$n.het.nb <- as.numeric(scale(etraits$n.het.nb))

# linear mixed effects model approach
efit <-lmer(log(emergence) ~ n.het.nb + n.con.nb + patch + elevation + (1|siteyear),data=etraits)
summary(efit)
res <- simulateResiduals(efit)
base::plot(res) # residuals bad

# quantile regression approach

# earliest emerging
e0.25 <- lqmm(
  fixed = emergence ~ n.het.nb + n.con.nb + patch + elevation,
  random = ~ 1,          # random intercept for site
  group = siteyear,
  tau = 0.25,             
  data = etraits
)
summary(e0.25)
# does this match quantreg?
# mod <- rq(emergence ~ n.het.nb + n.con.nb + patch + elevation,data=etraits,tau=0.25)
# summary(mod)
# Qualitatively yes, numbers close but not exact

# median
e0.5 <- lqmm(
  fixed = emergence ~ n.het.nb + n.con.nb + patch + elevation,
  random = ~ 1,          # random intercept for site
  group = siteyear,
  tau = 0.5,             
  data = etraits
)
summary(e0.5)
# latest emerging
e0.75 <- lqmm(
  fixed = emergence ~ n.het.nb + n.con.nb + patch + elevation,
  random = ~ 1,          # random intercept for site
  group = siteyear,
  tau = 0.75,             
  data = etraits
)
summary(e0.75)


## flowering phenology

ftraits <- data[-which(is.na(data$flower)),]
ftraits$flower <- as.numeric(strftime(ftraits$flower, format = "%V"))

# set up siteyear random effect & scaled predictors
ftraits$siteyear <- paste0(ftraits$site,ftraits$year)
ftraits$siteyear <- as.numeric(as.factor(ftraits$siteyear))
ftraits$n.con.nb <- as.numeric(scale(ftraits$n.con.nb))
ftraits$n.het.nb <- as.numeric(scale(ftraits$n.het.nb))

# linear mixed effects model appraoch
hist(ftraits$flower) # maybe normal? check residuals
ffit <-lmer(flower ~ n.het.nb + n.con.nb + patch + elevation + (1|siteyear),data=ftraits)
summary(ffit)
res <- simulateResiduals(ffit)
base::plot(res) # residuals bad

# quantile regression approach

# earliest flowering
f0.25 <- lqmm(
  fixed = flower ~ n.het.nb + n.con.nb + patch + elevation,
  random = ~ 1,          # random intercept for site
  group = siteyear,
  tau = 0.25,             
  data = ftraits
)
summary(f0.25)
# median
f0.5 <- lqmm(
  fixed = flower ~ n.het.nb + n.con.nb + patch + elevation,
  random = ~ 1,          # random intercept for site
  group = siteyear,
  tau = 0.5,             
  data = ftraits
)
summary(f0.5)
# latest flowering
f0.75 <- lqmm(
  fixed = flower ~ n.het.nb + n.con.nb + patch + elevation,
  random = ~ 1,          # random intercept for site
  group = siteyear,
  tau = 0.75,            
  data = ftraits
)
summary(f0.75)
