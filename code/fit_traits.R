### fitting models with traits as response

# Bayesian model for SLA, SRL, max height
trait.mod <- c("
data {
    int<lower=0> N; // number of observations
    real y[N]; // response variable
    int<lower=0> Nsiteyear; // number of site-year combinations
    int<lower=0,upper=Nsiteyear> siteyear[N]; // site-years 
    
    vector[N] e; // elevation covariate
    vector[N] p; // patch covariate
    vector[N] cn; // conspecific neighbor covariate
    vector[N] hn; // heterospecific neighbor covariate
}

parameters {
    vector[Nsiteyear] beta0; // hierarchical intercept
    vector[8] beta; // coefficients for covariates
    real alpha; // prior
    real<lower=0> nu; // prior
    real<lower=0> sigma; // standard deviation for sampling distribution
}    

transformed parameters {
  
  vector[N] mu; // storage of means
  
  for(i in 1:N){
  	  mu[i] = beta0[siteyear[i]] + beta[1]*e[i] + beta[2]*p[i] + beta[3]*cn[i] + beta[4]*hn[i] + beta[5]*e[i]*cn[i] + beta[6]*e[i]*hn[i] + beta[7]*p[i]*cn[i] + beta[8]*p[i]*hn[i];
  }
}

model {

	// Normal sampling distribution
  for(i in 1:N){
      y[i] ~ normal(mu[i],sigma);
  }
	
	// priors
	beta0 ~ normal(alpha,nu); // hierarchical intercept
	beta ~ normal(0,10);
	sigma ~ normal(0,10);
	
	// hyper priors
  alpha ~ normal(0,10); // alpha for beta0 prior
  nu ~ normal(0,10); // nu for beta0 prior
}
")

# set up stan
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# import data
data <- read.csv("./../clean_data/all_data_combined.csv")
data <- data[data$species=="BRTE",]

## SLA

# calculate SLA
ldata <- data[-which(is.na(data$total_leaf_area)),]
ldata$total_leaf_area <- ldata$total_leaf_area/100
ldata$SLA <- as.numeric(ldata$total_leaf_area)/as.numeric(ldata$total_leaf_mass)

# set up siteyear random effect & scaled predictors
ldata$siteyear <- paste0(ldata$site,ldata$year)
ldata$siteyear <- as.numeric(as.factor(ldata$siteyear))
ldata$elevation <- as.numeric(as.factor(ldata$elevation))
ldata$patch <- as.numeric(as.factor(ldata$patch))
ldata$n.con.nb <- as.numeric(scale(ldata$n.con.nb))
ldata$n.het.nb <- as.numeric(scale(ldata$n.het.nb))

# run model for SLA
sladata <- list(N=dim(ldata)[1],
             Nsiteyear=length(unique(ldata$siteyear)),
             siteyear=ldata$siteyear,
             y=ldata$SLA,
             e=ldata$elevation,
             p=ldata$patch,
             cn=ldata$n.con.nb,
             hn=ldata$n.het.nb)

SLA.fit <- stan(model_code = trait.mod,init=0,data=sladata)
summary(SLA.fit,pars=c("beta"))
plot(SLA.fit,pars=c("beta"))

# look at distribution of trait
# hist(log(ltraits$SLA)) # log-normal
# # model
# lfit <- lmer(log(SLA) ~ n.het.nb + n.con.nb + patch + elevation + (1|siteyear),data=ltraits)
# summary(lfit)
# res <- simulateResiduals(lfit)
base::plot(res) # residuals good enough
# positive effect of low elevation
# negative effect of heterospecifics
# positive effect of conspecifics


## SRL

# calculate SRL
rdata <- data[-which(is.na(data$length_cm)),]
rdata <- rdata[-which(is.na(rdata$root_mass)),]
rdata <- rdata[-which(rdata$root_mass==0),]
rdata$SRL <- as.numeric(rdata$length_cm)/as.numeric(rdata$root_mass)

# set up siteyear random effect & scaled predictors
rdata$siteyear <- paste0(rdata$site,rdata$year)
rdata$siteyear <- as.numeric(as.factor(rdata$siteyear))
rdata$elevation <- as.numeric(as.factor(rdata$elevation))
rdata$patch <- as.numeric(as.factor(rdata$patch))
rdata$n.con.nb <- as.numeric(scale(rdata$n.con.nb))
rdata$n.het.nb <- as.numeric(scale(rdata$n.het.nb))

# run model for SRL
srldata <- list(N=dim(rdata)[1],
                Nsiteyear=length(unique(rdata$siteyear)),
                siteyear=rdata$siteyear,
                y=rdata$SRL,
                e=rdata$elevation,
                p=rdata$patch,
                cn=rdata$n.con.nb,
                hn=rdata$n.het.nb)

SRL.fit <- stan(model_code = trait.mod,init=0,data=srldata)
summary(SRL.fit,pars=c("beta"))
plot(SRL.fit,pars=c("beta"))

# look at distribution of trait
# hist(log(rtraits$SRL)) # log-normal
# # fit model
# rfit <- lmer(log(SRL) ~ n.het.nb + n.con.nb + patch + elevation + (1|siteyear),data=rtraits)
# summary(rfit)
# res <- simulateResiduals(rfit)
# base::plot(res) # residuals good
# positive effect of conspecifics
# positive effect of shrub
# negative effect of low elevation


## height

hdata <- data[-which(is.na(data$max.height)),]
hdata <- hdata[-which(hdata$max.height==0),]

# set up siteyear random effect & scaled predictors
hdata$siteyear <- paste0(hdata$site,hdata$year)
hdata$siteyear <- as.numeric(as.factor(hdata$siteyear))
hdata$elevation <- as.numeric(as.factor(hdata$elevation))
hdata$patch <- as.numeric(as.factor(hdata$patch))
hdata$n.con.nb <- as.numeric(scale(hdata$n.con.nb))
hdata$n.het.nb <- as.numeric(scale(hdata$n.het.nb))

# run model for max height
maxhdata <- list(N=dim(hdata)[1],
                Nsiteyear=length(unique(hdata$siteyear)),
                siteyear=hdata$siteyear,
                y=hdata$max.height,
                e=hdata$elevation,
                p=hdata$patch,
                cn=hdata$n.con.nb,
                hn=hdata$n.het.nb)

height.fit <- stan(model_code = trait.mod,init=0,data=maxhdata)
summary(height.fit,pars=c("beta"))
plot(height.fit,pars=c("beta"))

# look at distribution of trait
# hist(log(htraits$max.height)) # log-normal
# # fit model
# hfit <-lmer(log(max.height) ~ n.het.nb + n.con.nb + patch + elevation + (1|siteyear),data=htraits)
# summary(hfit)
# res <- simulateResiduals(hfit)
# base::plot(res) # residuals good
# positive effect of conspecifics
# positive effect of shrubs
# psitivie effect of low elevation
