#### fitting models with biomass as response

# import data
data <- read.csv("./../clean_data/all_data_combined.csv")
data <- data[data$species=="BRTE",]

# formatting
bdata <- data[!is.na(data$aboveground_mass),]
bdata$siteyear <- paste0(bdata$site,bdata$year)
bdata$siteyear <- as.numeric(as.factor(bdata$siteyear))
bdata$nb.number <- bdata$n.con.nb + bdata$n.het.nb

bdata$elevation <- as.numeric(as.factor(bdata$elevation))
bdata$patch <- as.numeric(as.factor(bdata$patch))
# scaling
bdata$n.con.nb <- as.numeric(scale(bdata$n.con.nb))
bdata$n.het.nb <- as.numeric(scale(bdata$n.het.nb))

# Bayesian model

mass.mod <- c("
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

# adding data
data <- list(N=dim(bdata)[1],
             Nsiteyear=length(unique(bdata$siteyear)),
             siteyear=bdata$siteyear,
             y=bdata$aboveground_mass,
             e=bdata$elevation,
             p=bdata$patch,
             cn=bdata$n.con.nb,
             hn=bdata$n.het.nb)

mass.fit <- stan(model_code = mass.mod,init=0,data=data)
summary(mass.fit,pars=c("beta"))
plot(mass.fit,pars=c("beta"))

# model
mass <- lmer(log(aboveground_mass) ~ n.con.nb*patch + n.het.nb*patch + n.con.nb*elevation + n.het.nb*elevation + 
              (1|siteyear),data=bdata)
summary(mass)

# neighbor number model
# mass <- lmer(log(aboveground_mass) ~ elevation*nb.number + patch*nb.number+ 
#                (1|siteyear),data=biomass)
# summary(mass)

# quick plots
boxplot(mean.season.VWC ~ patch,data = biomass)
boxplot(mean.season.VWC ~ elevation,data = biomass)


