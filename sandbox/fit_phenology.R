## Fit emergence and flowering phenology data with linear and quantile regressions

pheno.mod <-c("
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

generated quantities{
    vector[N] res; // residuals
    vector[N] ypred; // replicated data

    for(i in 1:N){
    // sample replicated data
        ypred[i] = lognormal_rng(mu[i],sigma);
  
    // compute Pearson residuals
        res[i] = (y[i] - exp(mu[i]))/sqrt(exp(mu[i]));        
    }
}

")


# set up stan
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# import data
data <- read.csv("./../clean_data/all_data_combined.csv")
data <- data[data$species=="BRTE",]

## EMERGENCE PHENOLOGY

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
                  hn=etraits$n.het.nb)
# fit model
emerg.fit <- stan(model_code = pheno.mod,init=0,data=emod.data,iter=12000,warmup=6000)
save(emerg.fit,file=paste0("./../outputs/",Sys.Date(),"/stan_fits/emerg.fit.rda"))
# look at output
summary(emerg.fit,pars=c("beta"))
plot(emerg.fit,pars=c("beta"))

# save traceplot to show model convergence
trace.emerg <- traceplot(emerg.fit,pars=c("beta"))
trace.emerg
ggsave(trace.emerg,file=paste0(paste0("./../outputs/",Sys.Date(),"/trace_emerg.jpeg")),height = 6,width = 10)

## Posterior predictive checks

# get quantiles from posteriors
ypred<-extract(emerg.fit)$ypred
ypred_quant<-apply(ypred,2,quantile,probs=c(0.5,0.025,0.975))

# plot observed and predicted data
n<- length(etraits$emergence)

png(paste0("./../outputs/", Sys.Date(),"/emerg.ppc.png"), width = 8, height = 6, units = "in", res = 300)

plot(1:n,etraits$emergence,xlab="observation number",ylab="data value",pch=19)
segments(1:n,ypred_quant[2,],1:n,ypred_quant[3,],lty=3,col="firebrick")
points(1:n,ypred_quant[1,],pch=19,col=alpha("firebrick",.5))

dev.off()

# correlation between observed and predicted based on point estimates
cor(etraits$emergence,ypred_quant[1,])

# plot residuals
# png(paste0("./../outputs/", Sys.Date(),"/mass.resid.png"), width = 8, height = 6, units = "in", res = 300)
# 
# resid<-extract(mass.fit)$res
# plot(ypred_quant[1,],apply(resid,2,median),xlab="expected value",ylab="residual",pch=19)
# abline(h=0,lty=2)
# 
# dev.off()

### FLOWERING PHENOLOGY

ftraits <- data[-which(is.na(data$flower)),]
ftraits$flower <- as.numeric(strftime(ftraits$flower, format = "%V"))

# set up siteyear random effect & scaled predictors
ftraits$siteyear <- paste0(ftraits$year,ftraits$site,ftraits$shrub)
ftraits$siteyear <- as.numeric(as.factor(ftraits$siteyear))
ftraits$n.con.nb <- as.numeric(scale(ftraits$n.con.nb))
ftraits$n.het.nb <- as.numeric(scale(ftraits$n.het.nb))
ftraits$elevation <- as.numeric(as.factor(ftraits$elevation))
ftraits$patch <- as.numeric(as.factor(ftraits$patch))

# add data
emod.data <- list(N=dim(ftraits)[1],
                  Nsiteyear=length(unique(ftraits$siteyear)),
                  siteyear=ftraits$siteyear,
                  y=ftraits$flower,
                  e=ftraits$elevation,
                  p=ftraits$patch,
                  cn=ftraits$n.con.nb,
                  hn=ftraits$n.het.nb)
# fit model
flower.fit <- stan(model_code = pheno.mod,init=0,data=emod.data,iter=12000,warmup=6000)
save(flower.fit,file=paste0("./../outputs/",Sys.Date(),"/stan_fits/flower.fit.rda"))

# look at output
summary(flower.fit,pars=c("beta"))
plot(flower.fit,pars=c("beta"))

# save traceplot to show model convergence
trace.flower <- traceplot(flower.fit,pars=c("beta"))
trace.flower
ggsave(trace.flower,file=paste0(paste0("./../outputs/",Sys.Date(),"/trace_flower.jpeg")),height = 6,width = 10)

## Posterior predictive checks

# get quantiles from posteriors
ypred<-extract(flower.fit)$ypred
ypred_quant<-apply(ypred,2,quantile,probs=c(0.5,0.025,0.975))

# plot observed and predicted data
n<- length(ftraits$flower)

png(paste0("./../outputs/", Sys.Date(),"/flower.ppc.png"), width = 8, height = 6, units = "in", res = 300)

plot(1:n,ftraits$flower,xlab="observation number",ylab="data value",pch=19)
segments(1:n,ypred_quant[2,],1:n,ypred_quant[3,],lty=3,col="firebrick")
points(1:n,ypred_quant[1,],pch=19,col=alpha("firebrick",.5))

dev.off()

# correlation between observed and predicted based on point estimates
cor(ftraits$flower,ypred_quant[1,])

# plot residuals
# png(paste0("./../outputs/", Sys.Date(),"/mass.resid.png"), width = 8, height = 6, units = "in", res = 300)
# 
# resid<-extract(mass.fit)$res
# plot(ypred_quant[1,],apply(resid,2,median),xlab="expected value",ylab="residual",pch=19)
# abline(h=0,lty=2)
# 
# dev.off()

