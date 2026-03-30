#### fitting models with biomass as response

# import data
data <- read.csv("./../clean_data/all_data_combined.csv")
data <- data[data$species=="BRTE",]

# formatting
bdata <- data[!is.na(data$aboveground_mass),]
bdata$siteyear <- paste0(bdata$year,bdata$site,bdata$shrub)
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

# adding data
data <- list(N=dim(bdata)[1],
             Nsiteyear=length(unique(bdata$siteyear)),
             siteyear=bdata$siteyear,
             y=bdata$aboveground_mass,
             e=bdata$elevation,
             p=bdata$patch,
             cn=bdata$n.con.nb,
             hn=bdata$n.het.nb)

# run model
mass.fit <- stan(model_code = mass.mod,init=0,data=data,iter=12000,warmup=6000)

# look at output
summary(mass.fit,pars=c("beta"))
plot(mass.fit,pars=c("beta"))

# save traceplot to show model convergence
trace.mass <- traceplot(mass.fit,pars=c("beta"))
trace.mass
ggsave(trace.mass,file=paste0(paste0("./../outputs/",Sys.Date(),"/trace_mass.jpeg")),height = 6,width = 10)

## Posterior predictive checks

# get quantiles from posteriors
ypred<-extract(mass.fit)$ypred
ypred_quant<-apply(ypred,2,quantile,probs=c(0.5,0.025,0.975))

# plot observed and predicted data
n<- length(bdata$aboveground_mass)

png(paste0("./../outputs/", Sys.Date(),"/mass.ppc.png"), width = 8, height = 6, units = "in", res = 300)

plot(1:n,bdata$aboveground_mass,xlab="observation number",ylab="data value",pch=19)
segments(1:n,ypred_quant[2,],1:n,ypred_quant[3,],lty=3,col="firebrick")
points(1:n,ypred_quant[1,],pch=19,col=alpha("firebrick",.5))

dev.off()


## posterior predictive checks on SD
row_sds <- apply(ypred, 1, sd, na.rm = TRUE)
ypred.sd.lower <- quantile(row_sds,0.025) # bounds of 95% credible intervals
ypred.sd.upper <- quantile(row_sds,0.975)
true.sd <- sd(bdata$aboveground_mass)
# plot histogram of standard deviations of MCMC distributions
png(paste0("./../outputs/", Sys.Date(),"/mass.sd.ppc.png"), width = 8, height = 6, units = "in", res = 300)
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
true.mean <- mean(bdata$aboveground_mass)
# plot histogram of standard deviations of MCMC samples at each observation
png(paste0("./../outputs/", Sys.Date(),"/mass.mean.ppc.png"), width = 8, height = 6, units = "in", res = 300)
hist(row_means,xlab="Means of predicted y distribution",main="Posterior predictive checks of means")
abline(v=true.mean,col="red",lwd=2) # add true mean straight from data
abline(v=ypred.mean.lower,col="black",lty=3,lwd=2) # add 95% credible intervals
abline(v=ypred.mean.upper,col="black",lty=3,lwd=2)
# within 95% CIs
dev.off()

# correlation between observed and predicted based on point estimates
cor(bdata$aboveground_mass,ypred_quant[1,]) # 0.81

# plot residuals
png(paste0("./../outputs/", Sys.Date(),"/mass.resid.png"), width = 8, height = 6, units = "in", res = 300)

resid<-extract(mass.fit)$res
plot(ypred_quant[1,],apply(resid,2,median),xlab="expected value",ylab="residual",pch=19)
abline(h=0,lty=2)

dev.off()

t <- lmer(log(aboveground_mass) ~ elevation*n.con.nb + patch*n.con.nb + elevation*n.het.nb +
               patch*n.het.nb + (1|siteyear), data=bdata)
t <- lmer(log(aboveground_mass) ~ elevation*n.het.nb +
               patch*n.het.nb + (1|siteyear), data=bdata)
res <- resid(t)
plot(bdata$n.con.nb,res)

# save output
if(!dir.exists(paste0("./../outputs/",Sys.Date(),"/")))dir.create(paste0("./../outputs/",Sys.Date(),"/stan_fits/"))
save(mass.fit,file=paste0(paste0("./../outputs/",Sys.Date(),"/stan_fits/mass.fit.rda")))

