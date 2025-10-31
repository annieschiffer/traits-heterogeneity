#### fitting models with biomass as response

# import data
data <- read.csv("./../clean_data/all_data_combined.csv")
data <- data[data$species=="BRTE",]

biomass <- data[!is.na(data$aboveground_mass),]
biomass$siteyear <- paste0(biomass$site,biomass$year)
biomass$siteyear <- as.numeric(as.factor(biomass$siteyear))

# Bayesian model

biomass.mod <- c("
data {
    int<lower=0> N; // number of observations
    int<lower=1> y[N]; // response variable
    int<lower=0> Nsiteyear; // number of site-year combinations
    int<lower=0,upper=Nsiteyear> siteyear[N]; // site-years 
    
    vector[N] elev; // elevation covariate
    vector[N] ms; // microsite covariate
    vector[N] brte.nb; // number intraspecific neighbors covariate
    vector[N] inter.nb; // number interspecific neighbors covariate
}

parameters {
    vector[Nsiteyear] beta0; // hierarchical intercept
    vector[4] beta; // coefficients for covariates
    real<lower=0> sigma; // sigma
    real m; // prior
    real<lower=0> s; // prior
}    

model {
    
		// Normal sampling distribution and linear model
  	for(i in 1:N){
  	  mu[i] = beta0[siteyear[i]] + beta[1]*elev[i] + beta[2]*ms[i] + beta[3]*brte.nb[i] + beta[4]*inter.nb[i];
      y[i] ~ lognormal(mu[i],sigma); // lognormal data
      }
	
	// priors
	beta0 ~ normal(m,s); // hierarchical intercept
	beta ~ normal(0,10);
	sigma ~ normal(0,10);
	
	// hyper priors
  m ~ normal(0,10); // mu for beta0 prior
  s ~ normal(0,10); // sigma for beta0 prior
}
")

data <- list(y=biomass$aboveground_mass,
             N=length(biomass$aboveground_mass),
             Nsiteyear=length(unique(biomass$siteyear)),
             siteyear=biomass$siteyear,
             elev=biomass$elevation,
             ms=biomass$patch
             
  
)

