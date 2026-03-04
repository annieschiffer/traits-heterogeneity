### Fitting models with emergence as response

# import data
gdata <- read.csv("./../clean_data/all_data_combined.csv")
gdata <- gdata[gdata$species=="BRTE",]
gdata$nb.number <- gdata$n.con.nb + gdata$n.het.nb
# formatting
gdata$siteyear <- paste0(gdata$year,gdata$site,gdata$shrub)
gdata$siteyear <- as.numeric(as.factor(gdata$siteyear))

gdata$elevation <- as.numeric(as.factor(gdata$elevation))
gdata$patch <- as.numeric(as.factor(gdata$patch))

# scaling
gdata$nb.number <- as.numeric(scale(gdata$nb.number))
gdata$n.con.nb <- as.numeric(scale(gdata$n.con.nb))
gdata$n.het.nb <- as.numeric(scale(gdata$n.het.nb))

# Bayesian model

germ.mod <- c("
data {
    int<lower=0> N; // number of observations
    int<lower=0,upper=1> y[N]; // response variable
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
    real mu; // prior
    real<lower=0> sigma; // prior
}    

transformed parameters {
  
  vector[N] prob; // storage of probabilities
  
  for(i in 1:N){
  	  prob[i] = beta0[siteyear[i]] + beta[1]*e[i] + beta[2]*p[i] + beta[3]*cn[i] + beta[4]*hn[i] + beta[5]*e[i]*cn[i] + beta[6]*e[i]*hn[i] + beta[7]*p[i]*cn[i] + beta[8]*p[i]*hn[i];
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
             cn=gdata$n.con.nb,
             hn=gdata$n.het.nb)

# run model
germ.fit <- stan(model_code = germ.mod,init=0,data=data,iter=12000,warmup=6000)

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

# correlation between observed and predicted based on point estimates
cor(gdata$germination,ypred_quant[1,])

# plot residuals
png(paste0("./../outputs/", Sys.Date(),"/germ.resid.png"), width = 8, height = 6, units = "in", res = 300)

resid<-extract(germ.fit)$res
plot(ypred_quant[1,],apply(resid,2,median),xlab="expected value",ylab="residual",pch=19)
abline(h=0,lty=2)

dev.off()


# save output
save(germ.fit,file=paste0(paste0("./../outputs/",Sys.Date(),"/stan_fits/germ.fit.rda")))

