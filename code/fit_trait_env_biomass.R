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

#write.csv(bdata.complete,"./../clean_data/complete_data_biomass.csv",row.names = FALSE)

# scaling
bdata.complete$n.con.nb <- as.numeric(scale(bdata.complete$n.con.nb))
bdata.complete$n.het.nb <- as.numeric(scale(bdata.complete$n.het.nb))
bdata.complete$max.height <- as.numeric(scale(bdata.complete$max.height))
bdata.complete$SLA <- as.numeric(scale(bdata.complete$SLA))
bdata.complete$SRL <- as.numeric(scale(bdata.complete$SRL))

# formatting for stan
bdata.complete$elevation <- as.numeric(as.factor(bdata.complete$elevation))
bdata.complete$patch <- as.numeric(as.factor(bdata.complete$patch))

bdata.complete$siteyear <- paste0(bdata.complete$year,bdata.complete$site,bdata.complete$shrub)
bdata.complete$siteyear <- as.numeric(as.factor(bdata.complete$siteyear))

# Bayesian model

TbyE.mass.mod <- c("
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
      y[i] ~ lognormal(mu[i],sigma);
  }
	
	// priors
	beta0 ~ normal(alpha,nu); // hierarchical intercept
	beta ~ normal(0,5);
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

# run model
TbyE.mass.fit <- stan(model_code = TbyE.mass.mod,init=0,data=data,iter = 12000,warmup = 6000)

# look at output
summary(TbyE.mass.fit,pars=c("beta"))
plot(TbyE.mass.fit,pars=c("beta"))

# save traceplot to show model convergence
trace.TbyE.mass <- traceplot(TbyE.mass.fit,pars=c("beta"))
trace.TbyE.mass
ggsave(trace.TbyE.mass,file=paste0(paste0("./../outputs/",Sys.Date(),"/trace_TbyE.mass.jpeg")),height = 6,width = 10)

## Posterior predictive checks

# get quantiles from posteriors
ypred<-extract(TbyE.mass.fit)$ypred
ypred_quant<-apply(ypred,2,quantile,probs=c(0.5,0.025,0.975))

# plot observed and predicted data
n<- length(bdata.complete$aboveground_mass)

png(paste0("./../outputs/", Sys.Date(),"/TbyE.mass.ppc.png"), width = 8, height = 6, units = "in", res = 300)

plot(1:n,bdata.complete$aboveground_mass,xlab="observation number",ylab="data value",pch=19)
segments(1:n,ypred_quant[2,],1:n,ypred_quant[3,],lty=3,col="firebrick")
points(1:n,ypred_quant[1,],pch=19,col=alpha("firebrick",.5))

dev.off()

## posterior predictive checks on SD
row_sds <- apply(ypred, 1, sd, na.rm = TRUE)
ypred.sd.lower <- quantile(row_sds,0.025) # bounds of 95% credible intervals
ypred.sd.upper <- quantile(row_sds,0.975)
true.sd <- sd(bdata.complete$aboveground_mass)
# plot histogram of standard deviations of MCMC distributions
png(paste0("./../outputs/", Sys.Date(),"/TbyE.sd.ppc.png"), width = 8, height = 6, units = "in", res = 300)
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
true.mean <- mean(bdata.complete$aboveground_mass)
# plot histogram of standard deviations of MCMC samples at each observation
png(paste0("./../outputs/", Sys.Date(),"/TbyE.mean.ppc.png"), width = 8, height = 6, units = "in", res = 300)
hist(row_means,xlab="Means of predicted y distribution",main="Posterior predictive checks of means")
abline(v=true.mean,col="red",lwd=2) # add true mean straight from data
abline(v=ypred.mean.lower,col="black",lty=3,lwd=2) # add 95% credible intervals
abline(v=ypred.mean.upper,col="black",lty=3,lwd=2)
# within 95% CIs
dev.off()

# correlation between observed and predicted based on point estimates
cor(bdata.complete$aboveground_mass,ypred_quant[1,]) # 0.81

# plot residuals
png(paste0("./../outputs/", Sys.Date(),"/TbyE.mass.resid.png"), width = 8, height = 6, units = "in", res = 300)

resid<-extract(TbyE.mass.fit)$res
plot(ypred_quant[1,],apply(resid,2,median),xlab="expected value",ylab="residual",pch=19)
abline(h=0,lty=2)

dev.off()

# save output
if(!dir.exists(paste0("./../outputs/",Sys.Date(),"/")))dir.create(paste0("./../outputs/",Sys.Date(),"/stan_fits/"))
save(TbyE.mass.fit,file=paste0(paste0("./../outputs/",Sys.Date(),"/stan_fits/TbyE.mass.fit.rda")))
