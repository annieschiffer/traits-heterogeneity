### Fitting models with emergence as response

# import data
gdata <- read.csv("./../clean_data/all_data_combined.csv")
gdata <- gdata[gdata$species=="BRTE",]

# formatting
gdata$siteyear <- paste0(gdata$year,gdata$site,gdata$shrub)
gdata$siteyear <- as.numeric(as.factor(gdata$siteyear))

gdata$elevation <- as.numeric(as.factor(gdata$elevation))
gdata$patch <- as.numeric(as.factor(gdata$patch))
gdata$treat <- as.numeric(as.factor(gdata$subplot))

# scaling
# gdata$nb.number <- gdata$n.con.nb + gdata$n.het.nb
# gdata$nb.number <- as.numeric(scale(gdata$nb.number))
# gdata$n.con.nb <- as.numeric(scale(gdata$n.con.nb))
# gdata$n.het.nb <- as.numeric(scale(gdata$n.het.nb))

# Bayesian model

germ.mod <- c("
data {
    int<lower=0> N; // number of observations
    int<lower=0,upper=1> y[N]; // response variable
    int<lower=0> Nsiteyear; // number of site-year combinations
    int<lower=0,upper=Nsiteyear> siteyear[N]; // site-years 
    
    vector[N] e; // elevation covariate
    vector[N] p; // patch covariate
    vector[N] t; // competition treatment covariate
}

parameters {
    vector[Nsiteyear] beta0; // hierarchical intercept
    vector[5] beta; // coefficients for covariates
    real mu; // prior
    real<lower=0> sigma; // prior
}    

transformed parameters {
  
  vector[N] prob; // storage of probabilities
  
  for(i in 1:N){
  	  prob[i] = beta0[siteyear[i]] + beta[1]*e[i] + beta[2]*p[i] + beta[3]*t[i] + beta[4]*e[i]*t[i] + beta[5]*p[i]*t[i];
  }
}

model {

	// Bernoulli sampling distribution and linear model
  for(i in 1:N){
      y[i] ~ bernoulli_logit(prob[i]);
  }
	
	// priors
	beta0 ~ normal(mu,sigma); // hierarchical intercept
	beta ~ normal(0,10);
	
	// hyper priors
  mu ~ normal(0,10); // mu for beta0 prior
  sigma ~ normal(0,10); // sigma for beta0 prior
}

generated quantities{
    vector[N] res; // residuals
    vector[N] ypred; // replicated data

    for(i in 1:N){
    // sample replicated data
        ypred[i] = bernoulli_rng(inv_logit(prob[i]));
  
    // compute Pearson residuals
        res[i] = (y[i] - inv_logit(prob[i]))/sqrt(inv_logit(prob[i]));        
    }
}

")

# set up stan
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# adding data
data <- list(N=dim(gdata)[1],
             Nsiteyear=length(unique(gdata$siteyear)),
             siteyear=gdata$siteyear,
             y=gdata$germination,
             e=gdata$elevation,
             p=gdata$patch,
             t=gdata$treat)

# run model
germ.fit <- stan(model_code = germ.mod,init=0,data=data, iter=12000, warmup=6000)

# look at output
summary(germ.fit,pars=c("beta"))
plot(germ.fit,pars=c("beta"))

# save traceplot to show model convergence
trace.germ <- traceplot(germ.fit,pars=c("beta"))
trace.germ
ggsave(trace.germ,file=paste0(paste0("./../outputs/",Sys.Date(),"/trace_germ.jpeg")),height = 6,width = 10)


## Posterior predictive checks

# get quantiles from posteriors
ypred<-extract(germ.fit)$ypred
ypred_quant<-apply(ypred,2,quantile,probs=c(0.5,0.025,0.975))

# plot observed and predicted data
n<- length(gdata$germination)

png(paste0("./../outputs/", Sys.Date(),"/germ.ppc.png"), width = 8, height = 6, units = "in", res = 300)

plot(1:n,gdata$germination,xlab="observation number",ylab="data value",pch=19)
#segments(1:n,ypred_quant[2,],1:n,ypred_quant[3,],lty=3,col="firebrick")
points(1:n,ypred_quant[1,],pch=19,col=alpha("firebrick",.5))

dev.off()

## posterior predictive checks on SD
row_sds <- apply(ypred, 1, sd, na.rm = TRUE)
ypred.sd.lower <- quantile(row_sds,0.025) # bounds of 95% credible intervals
ypred.sd.upper <- quantile(row_sds,0.975)
true.sd <- sd(gdata$germination)
# plot histogram of standard deviations of MCMC distributions
png(paste0("./../outputs/", Sys.Date(),"/germ.sd.ppc.png"), width = 8, height = 6, units = "in", res = 300)
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
true.mean <- mean(gdata$germination)
# plot histogram of standard deviations of MCMC samples at each observation
png(paste0("./../outputs/", Sys.Date(),"/germ.mean.ppc.png"), width = 8, height = 6, units = "in", res = 300)
hist(row_means,xlab="Means of predicted y distribution",main="Posterior predictive checks of means")
abline(v=true.mean,col="red",lwd=2) # add true mean straight from data
abline(v=ypred.mean.lower,col="black",lty=3,lwd=2) # add 95% credible intervals
abline(v=ypred.mean.upper,col="black",lty=3,lwd=2)
# within 95% CIs
dev.off()

# correlation between observed and predicted based on point estimates
cor(gdata$germination,ypred_quant[1,])

# plot residuals
png(paste0("./../outputs/", Sys.Date(),"/germ.resid.png"), width = 8, height = 6, units = "in", res = 300)

resid<-extract(germ.fit)$res
plot(ypred_quant[1,],apply(resid,2,median),xlab="expected value",ylab="residual",pch=19)
abline(h=0,lty=2)

dev.off()
# 
# t <- glmmTMB(germination ~ elevation*n.con.nb + patch*n.con.nb + elevation*n.het.nb +
#                patch*n.het.nb + (1|siteyear), data=gdata, family = binomial(link="logit"))
# t <- glmmTMB(germination ~ elevation*n.het.nb +
#                patch*n.het.nb + (1|siteyear), data=gdata, family = binomial(link="logit"))

# save output
save(germ.fit,file=paste0(paste0("./../outputs/",Sys.Date(),"/stan_fits/germ.fit.rda")))

