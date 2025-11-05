#### fitting models with biomass as response

# import data
data <- read.csv("./../clean_data/all_data_combined.csv")
data <- data[data$species=="BRTE",]

biomass <- data[!is.na(data$aboveground_mass),]
biomass$siteyear <- paste0(biomass$site,biomass$year)
biomass$siteyear <- as.numeric(as.factor(biomass$siteyear))
biomass$elevation <- ifelse(biomass$elevation == "low", 0, 1)
biomass$patch <- ifelse(biomass$patch == "open", 0, 1)

# Bayesian model

biomass.mod <- c("
data {
    int<lower=0> N; // number of observations
    real y[N]; // response variable
    int<lower=0> Nsiteyear; // number of site-year combinations
    int<lower=0,upper=Nsiteyear> siteyear[N]; // site-years 
    
    vector[N] elev; // elevation covariate
    vector[N] ms; // microsite covariate
    vector[N] con; // number conspecific neighbors covariate
    vector[N] het; // number heterospecific neighbors covariate
}

parameters {
    vector[Nsiteyear] beta0; // hierarchical intercept
    vector[4] beta; // coefficients for covariates
    real<lower=0> sigma; // sigma
    real m; // prior
    real<lower=0> s; // prior
}    

model {

    real mu; // mean
    
		// Normal sampling distribution and linear model
  	for(i in 1:N){
  	  mu = beta0[siteyear[i]] + beta[1]*elev[i] + beta[2]*ms[i] + beta[3]*con[i] + beta[4]*het[i];
      y[i] ~ lognormal(mu,sigma); // lognormal data
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

# set up stan
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

data <- list(y=biomass$aboveground_mass,
             N=length(biomass$aboveground_mass),
             Nsiteyear=length(unique(biomass$siteyear)),
             siteyear=biomass$siteyear,
             elev=biomass$elevation,
             ms=biomass$patch,
             con=biomass$n.con.nb,
             het=biomass$n.het.nb)

biomass.fit <- stan(model_code = biomass.mod,data=data,iter=12000,warmup=6000)
biomass.fit

lmer(log(aboveground_mass) ~ n.con.nb + n.het.nb + patch + elevation + (1|siteyear),data=biomass)
# exactly the same as bayesian model
