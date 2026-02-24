### Fitting models with aboveground biomass as a function of trait x environment

# import data
data <- read.csv("./../clean_data/all_data_combined.csv")
data <- data[data$species=="BRTE",]

# only complete data set
bdata <- data[!is.na(data$aboveground_mass),]
traits <- bdata[,c("max.height","length_cm","root_mass","total_leaf_area","total_leaf_mass")]
bdata.complete <- bdata[which(complete.cases(traits)==TRUE),]

# calculate SLA
bdata.complete$total_leaf_area <- bdata.complete$total_leaf_area/100
bdata.complete$SLA <- as.numeric(bdata.complete$total_leaf_area)/as.numeric(bdata.complete$total_leaf_mass)

# calculate SRL
bdata.complete$SRL <- as.numeric(bdata.complete$length_cm)/as.numeric(bdata.complete$root_mass)
bdata.complete <- bdata.complete[-which(bdata.complete$SRL==Inf),]

# scaling
bdata.complete$n.con.nb <- as.numeric(scale(bdata.complete$n.con.nb))
bdata.complete$n.het.nb <- as.numeric(scale(bdata.complete$n.het.nb))
bdata.complete$max.height <- as.numeric(scale(bdata.complete$max.height))
bdata.complete$SLA <- as.numeric(scale(bdata.complete$SLA))
bdata.complete$SRL <- as.numeric(scale(bdata.complete$SRL))

# formatting for stan
#bdata.complete$env <- paste0(bdata.complete$elevation,"-",bdata.complete$patch)
bdata.complete$elevation <- as.numeric(as.factor(bdata.complete$elevation))
bdata.complete$patch <- as.numeric(as.factor(bdata.complete$patch))

bdata.complete$siteyear <- paste0(bdata.complete$site,bdata.complete$year)
bdata.complete$siteyear <- as.numeric(as.factor(bdata.complete$siteyear))

# Bayesian model

TbyE.mod <- c("
data {
    int<lower=0> N; // number of observations
    real y[N]; // response variable
    int<lower=0> Nsiteyear; // number of site-year combinations
    int<lower=0,upper=Nsiteyear> siteyear[N]; // site-years 
    
    vector[N] e; // elevation covariate
    vector[N] p; // patch covariate
    vector[N] SRL; // SRL covariate
    vector[N] SLA; // SLA covariate
    vector[N] h; // height covariate
}

parameters {
    vector[Nsiteyear] beta0; // hierarchical intercept
    vector[11] beta; // coefficients for covariates
    real alpha; // prior
    real<lower=0> nu; // prior
    real<lower=0> sigma; // standard deviation for sampling distribution
}    

transformed parameters {
  
  vector[N] mu; // storage of means
  
  for(i in 1:N){
  	  mu[i] = beta0[siteyear[i]] + beta[1]*e[i] + beta[2]*p[i] + beta[3]*SRL[i] + beta[4]*SLA[i] + beta[5]*h[i] + beta[6]*e[i]*SRL[i] + beta[7]*e[i]*SLA[i] + beta[8]*e[i]*h[i] + beta[9]*p[i]*SRL[i] + beta[10]*p[i]*SLA[i] + beta[11]*p[i]*h[i];
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

# adding data
data <- list(N=dim(bdata.complete)[1],
             Nsiteyear=length(unique(bdata.complete$siteyear)),
             siteyear=bdata.complete$siteyear,
             y=bdata.complete$aboveground_mass,
             SLA=bdata.complete$SLA,
             SRL=bdata.complete$SRL,
             h=bdata.complete$max.height,
             e=bdata.complete$elevation,
             p=bdata.complete$patch
)

TbyE.fit <- stan(model_code = TbyE.mod,init=0,data=data,iter = 4000,warmup = 2000)
summary(TbyE.fit,pars=c("beta"))
plot(TbyE.fit,pars=c("beta"))

# biomass as a function of traits, environmental conditions, and their interaction
biomass.trait <- lmer(log(aboveground_mass) ~ env*SLA + env*SRL + env*max.height + env*emergence +
                (1|siteyear),data=complete)
summary(biomass.trait)
