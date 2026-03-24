## Fit emergence and flowering phenology data with linear and quantile regressions

### quantile regression or gamma?? 

pheno.mod <-c("

data {
    int<lower=0> N; // number of observations
    real y[N]; // response variable
    int<lower=0> Nsiteyear; // number of site-year combinations
    int<lower=0,upper=Nsiteyear> siteyear[N]; // site-years 
    real<lower=0, upper=1> tau; // quantile
    
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
      y[i] ~ skew_double_exponential(mu[i],sigma,tau);
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

## emergence phenology

etraits <- data[-which(is.na(data$emergence)),]
etraits$emergence <- as.numeric(strftime(etraits$emergence, format = "%V"))

# set up siteyear random effect & scaled predictors
etraits$siteyear <- paste0(etraits$year,etraits$site,etraits$shrub)
etraits$siteyear <- as.numeric(as.factor(etraits$siteyear))
etraits$n.con.nb <- as.numeric(scale(etraits$n.con.nb))
etraits$n.het.nb <- as.numeric(scale(etraits$n.het.nb))
etraits$elevation <- as.numeric(as.factor(etraits$elevation))
etraits$patch <- as.numeric(as.factor(etraits$patch))

# add data
emod.data <- list(N=dim(etraits)[1],
                  Nsiteyear=length(unique(etraits$siteyear)),
                  siteyear=etraits$siteyear,
                  y=etraits$emergence,
                  e=etraits$elevation,
                  p=etraits$patch,
                  cn=etraits$n.con.nb,
                  hn=etraits$n.het.nb,
                  tau=0.25)
# fit model
e25.fit <- stan(model_code = pheno.mod,init=0,data=emod.data,iter=12000,warmup=6000)
summary(e25.fit,pars=c("beta"))
plot(e25.fit,pars=c("beta"))



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
