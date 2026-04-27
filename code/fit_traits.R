### fitting models with traits as response

# Bayesian model for SLA, SRL, max height, emergence and flowering day of year
trait.mod <- c("
data {
    int<lower=0> N; // number of observations
    real y[N]; // response variable
    int<lower=0> Nsiteyear; // number of site-year combinations
    int<lower=0,upper=Nsiteyear> siteyear[N]; // site-years 
    
    vector[N] e; // elevation covariate
    vector[N] p; // patch covariate
    vector[N] t; // competition treatment covariate
}

parameters {
    vector[Nsiteyear] beta0; // hierarchical intercept
    vector[5] beta; // coefficients for covariates
    real alpha; // prior
    real<lower=0> nu; // prior
    real<lower=0> sigma; // standard deviation for sampling distribution
}    

transformed parameters {
  
  vector[N] mu; // storage of means
  
  for(i in 1:N){
  	  mu[i] = beta0[siteyear[i]] + beta[1]*e[i] + beta[2]*p[i] + beta[3]*t[i] + beta[4]*e[i]*t[i] + beta[5]*p[i]*t[i];
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
  if(trait=="emerg"){
    tdata <- all.data[-which(is.na(all.data$emergence)),]
    tdata$focal.trait <- as.numeric(strftime(tdata$emergence, format = "%V"))
  }
  if(trait=="flower"){
    tdata <- all.data[-which(is.na(all.data$flower)),]
    tdata$focal.trait <- as.numeric(strftime(tdata$flower, format = "%V"))
  }
  # set up siteyear random effect & scaled predictors
  tdata$siteyear <- paste0(tdata$year,tdata$site,tdata$shrub)
  tdata$siteyear <- as.numeric(as.factor(tdata$siteyear))
  tdata$elevation <- as.numeric(as.factor(tdata$elevation)) 
  tdata$patch <- as.numeric(as.factor(tdata$patch))
  tdata$treat <- as.numeric(as.factor(tdata$subplot))
  
  # run model for SLA
  tmod.data <- list(N=dim(tdata)[1],
                  Nsiteyear=length(unique(tdata$siteyear)),
                  siteyear=tdata$siteyear,
                  y=tdata$focal.trait,
                  e=tdata$elevation,
                  p=tdata$patch,
                  t=tdata$treat)
  return(tmod.data)
}

# fit SLA model
sladata <- format.trait("SLA",data)
SLA.fit <- stan(model_code = trait.mod,init=0,data=sladata,iter=12000,warmup=6000)
summary(SLA.fit,pars=c("beta"))
plot(SLA.fit,pars=c("beta"))
# save traceplot to show model convergence
trace.SLA <- traceplot(SLA.fit,pars=c("beta"))
trace.SLA
ggsave(trace.SLA,file=paste0(paste0("./../outputs/",Sys.Date(),"/trace_SLA.jpeg")),height = 6,width = 10)

# fit SRL model
srldata <- format.trait("SRL",data)
SRL.fit <- stan(model_code = trait.mod,init=0,data=srldata,iter=12000,warmup=6000)
summary(SRL.fit,pars=c("beta"))
plot(SRL.fit,pars=c("beta"))
# save traceplot to show model convergence
trace.SRL <- traceplot(SRL.fit,pars=c("beta"))
trace.SRL
ggsave(trace.SRL,file=paste0(paste0("./../outputs/",Sys.Date(),"/trace_SRL.jpeg")),height = 6,width = 10)
# when competitors present, 

# fit height model 
maxhdata <- format.trait("max.height",data)
height.fit <- stan(model_code = trait.mod,init=0,data=maxhdata,iter=12000,warmup=6000)
summary(height.fit,pars=c("beta"))
plot(height.fit,pars=c("beta"))
# save traceplot to show model convergence
trace.height <- traceplot(height.fit,pars=c("beta"))
trace.height
ggsave(trace.height,file=paste0(paste0("./../outputs/",Sys.Date(),"/trace_height.jpeg")),height = 6,width = 10)

# fit emergence phenology model 
edata <- format.trait("emerg",data)
emerg.fit <- stan(model_code = trait.mod,init=0,data=edata,iter=12000,warmup=6000)
summary(emerg.fit,pars=c("beta"))
plot(emerg.fit,pars=c("beta"))
# save traceplot to show model convergence
trace.emerg <- traceplot(emerg.fit,pars=c("beta"))
trace.emerg
ggsave(trace.emerg,file=paste0(paste0("./../outputs/",Sys.Date(),"/trace_emerg.jpeg")),height = 6,width = 10)

# fit flowering phenology model 
fdata <- format.trait("flower",data)
flower.fit <- stan(model_code = trait.mod,init=0,data=fdata,iter=12000,warmup=6000)
summary(flower.fit,pars=c("beta"))
plot(flower.fit,pars=c("beta"))
# save traceplot to show model convergence
trace.flower <- traceplot(flower.fit,pars=c("beta"))
trace.flower
ggsave(trace.flower,file=paste0(paste0("./../outputs/",Sys.Date(),"/trace_flower.jpeg")),height = 6,width = 10)

## Posterior predictive checks

get.ppc <- function(stanfit,trait,response){
  
  # get quantiles from posteriors
  ypred<-extract(stanfit)$ypred
  ypred_quant<-apply(ypred,2,quantile,probs=c(0.5,0.025,0.975))
  
  # plot observed and predicted data
  n<- length(response)
  
  png(paste0("./../outputs/", Sys.Date(),"/",trait,".ppc.png"), width = 8, height = 6, units = "in", res = 300)
  
  plot(1:n,response,xlab="observation number",ylab="data value",pch=19)
  segments(1:n,ypred_quant[2,],1:n,ypred_quant[3,],lty=3,col="firebrick")
  points(1:n,ypred_quant[1,],pch=19,col=alpha("firebrick",.5))
  dev.off()
  
  ## posterior predictive checks on SD
  row_sds <- apply(ypred, 1, sd, na.rm = TRUE)
  ypred.sd.lower <- quantile(row_sds,0.025) # bounds of 95% credible intervals
  ypred.sd.upper <- quantile(row_sds,0.975)
  true.sd <- sd(response)
  # plot histogram of standard deviations of MCMC distributions
  png(paste0("./../outputs/", Sys.Date(),"/",trait,".sd.ppc.png"), width = 8, height = 6, units = "in", res = 300)
  hist(row_sds,xlab="Standard deviation of predicted y distributions",main="Posterior predictive checks of standard deviations")
  abline(v=true.sd,col="red",lwd=2) # add true sd straight from data
  abline(v=ypred.sd.lower,col="black",lty=3,lwd=2) # add 95% credible intervals
  abline(v=ypred.sd.upper,col="black",lty=3,lwd=2)
  # within 95% CIs
  dev.off()
  
  ## posterior predictive checks on means
  row_means <- apply(ypred,1,mean,na.rm=TRUE)
  ypred.mean.lower <- quantile(row_means,0.025) # bounds of 95% credible intervals
  ypred.mean.upper <- quantile(row_means,0.975)
  true.mean <- mean(response)
  # plot histogram of standard deviations of MCMC samples at each observation
  png(paste0("./../outputs/", Sys.Date(),"/",trait,".mean.ppc.png"), width = 8, height = 6, units = "in", res = 300)
  hist(row_means,xlab="Means of predicted y distribution",main="Posterior predictive checks of means")
  abline(v=true.mean,col="red",lwd=2) # add true mean straight from data
  abline(v=ypred.mean.lower,col="black",lty=3,lwd=2) # add 95% credible intervals
  abline(v=ypred.mean.upper,col="black",lty=3,lwd=2)
  # within 95% CIs
  dev.off()
  
}

# load all models
# load("./../outputs/2026-03-04/stan_fits/SLA.fit.rda")
# load("./../outputs/2026-03-04/stan_fits/SRL.fit.rda")
# load("./../outputs/2026-03-04/stan_fits/height.fit.rda")
# load("./../outputs/2026-03-24/stan_fits/emerg.fit.rda")
# load("./../outputs/2026-03-24/stan_fits/flower.fit.rda")

get.ppc(SLA.fit,"SLA",sladata[["y"]])
get.ppc(SRL.fit,"SRL",srldata[["y"]])
get.ppc(height.fit,"height",maxhdata[["y"]])
get.ppc(emerg.fit,"emerg",edata[["y"]])
get.ppc(flower.fit,"flower",fdata[["y"]])


# save output
if(!dir.exists(paste0("./../outputs/",Sys.Date(),"/")))dir.create(paste0("./../outputs/",Sys.Date(),"/stan_fits/"))
save(SLA.fit,file=paste0(paste0("./../outputs/",Sys.Date(),"/stan_fits/SLA.fit.rda")))
save(SRL.fit,file=paste0(paste0("./../outputs/",Sys.Date(),"/stan_fits/SRL.fit.rda")))
save(height.fit,file=paste0(paste0("./../outputs/",Sys.Date(),"/stan_fits/height.fit.rda")))
save(emerg.fit,file=paste0(paste0("./../outputs/",Sys.Date(),"/stan_fits/emerg.fit.rda")))
save(flower.fit,file=paste0(paste0("./../outputs/",Sys.Date(),"/stan_fits/flower.fit.rda")))
