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
      y[i] ~ lognormal(mu[i],sigma);
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

# create function
format.trait <- function(trait,all.data){
  if(trait=="SLA"){
    # calculate SLA
    tdata <- all.data[-which(is.na(all.data$total_leaf_area)),]
    tdata$total_leaf_area <- tdata$total_leaf_area/100
    tdata$focal.trait <- as.numeric(tdata$total_leaf_area)/as.numeric(tdata$total_leaf_mass)
  }
  if(trait=="SRL"){
    # calculate SRL
    tdata <- all.data[-which(is.na(all.data$length_cm)),]
    tdata <- tdata[-which(is.na(tdata$root_mass)),]
    tdata <- tdata[-which(tdata$root_mass==0),]
    tdata$focal.trait <- as.numeric(tdata$length_cm)/as.numeric(tdata$root_mass)
  }
  if(trait=="max.height"){
    tdata <- all.data[-which(is.na(all.data$max.height)),]
    tdata <- tdata[-which(tdata$max.height==0),]
    tdata$focal.trait <- tdata$max.height
  }
  # set up siteyear random effect & scaled predictors
  tdata$siteyear <- paste0(tdata$year,tdata$site,tdata$shrub)
  tdata$siteyear <- as.numeric(as.factor(tdata$siteyear))
  tdata$elevation <- as.numeric(as.factor(tdata$elevation))
  tdata$patch <- as.numeric(as.factor(tdata$patch))
  tdata$n.con.nb <- as.numeric(scale(tdata$n.con.nb))
  tdata$n.het.nb <- as.numeric(scale(tdata$n.het.nb))
  
  # run model for SLA
  tmod.data <- list(N=dim(tdata)[1],
                  Nsiteyear=length(unique(tdata$siteyear)),
                  siteyear=tdata$siteyear,
                  y=tdata$focal.trait,
                  e=tdata$elevation,
                  p=tdata$patch,
                  cn=tdata$n.con.nb,
                  hn=tdata$n.het.nb)
  return(tmod.data)
}

# fit SLA model
sladata <- format.trait("SLA",data)
SLA.fit <- stan(model_code = trait.mod,init=0,data=sladata,iter=12000,warmup=6000)
summary(SLA.fit,pars=c("beta"))
plot(SLA.fit,pars=c("beta"))

# fit SRL model
srldata <- format.trait("SRL",data)
SRL.fit <- stan(model_code = trait.mod,init=0,data=srldata,iter=12000,warmup=6000)
summary(SRL.fit,pars=c("beta"))
plot(SRL.fit,pars=c("beta"))

# fit height model 
maxhdata <- format.trait("max.height",data)
height.fit <- stan(model_code = trait.mod,init=0,data=maxhdata,iter=12000,warmup=6000)
summary(height.fit,pars=c("beta"))
plot(height.fit,pars=c("beta"))

# save output
if(!dir.exists(paste0("./../outputs/",Sys.Date(),"/")))dir.create(paste0("./../outputs/",Sys.Date(),"/stan_fits/"))
save(SLA.fit,file=paste0(paste0("./../outputs/",Sys.Date(),"/stan_fits/SLA.fit.rda")))
save(SRL.fit,file=paste0(paste0("./../outputs/",Sys.Date(),"/stan_fits/SRL.fit.rda")))
save(height.fit,file=paste0(paste0("./../outputs/",Sys.Date(),"/stan_fits/height.fit.rda")))

