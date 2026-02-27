### Writing custom functions to format figure data

predict.emerg <- function(stanfit,emergdata){

  # extract the posterior samples for all betas
  post <- extract(stanfit)
  beta_elev <- post$beta[,1]
  beta_p <- post$beta[,2]
  beta_con <- post$beta[,3]
  beta_het <- post$beta[,4]
  beta_elxcon <- post$beta[,5]
  beta_elxhet <- post$beta[,6]
  beta_pxcon <- post$beta[,7]
  beta_pxhet <- post$beta[,8]
  
  # find the mean and standard deviation of raw data (I scaled data in stan model)
  mean.raw.con <- mean(emergdata$n.con.nb)
  mean.raw.het <- mean(emergdata$n.het.nb)
  sd.raw.con <- sd(emergdata$n.con.nb)
  sd.raw.het <- sd(emergdata$n.het.nb)
  
  # create new y data frame
  newdata <- expand.grid(
    con.raw = 0:50,
    het.raw = 0:50,
    elevation=c(0,1),
    patch=c(0,1)
  )
  
  # scale the number of neighbors, like I did for stan model
  newdata$con.scale <- (newdata$con.raw - mean.raw.con) / sd.raw.con
  newdata$het.scale <- (newdata$het.raw - mean.raw.het) / sd.raw.het
  
  # set up prediction matrix & prepare for loop
  inv_logit <- function(x){1/(1+exp(-x))}
  ndraws <- length(post$mu)
  npoints <- nrow(newdata)
  p_mat <- matrix(NA,ndraws,npoints)
  
  # pull beta posterior samples and use new data to estimate over new neighbor abundances
  for(i in 1:ndraws){
    eta <- post$mu[i] + beta_elev[i]*newdata$elevation + beta_p[i]*newdata$patch + beta_con[i]*newdata$con.scale + beta_het[i]*newdata$het.scale +
      beta_elxcon[i]*newdata$con.scale*newdata$elevation + beta_elxhet[i]*newdata$het.scale*newdata$elevation+
      beta_pxcon[i]*newdata$con.scale*newdata$patch + beta_pxhet[i]*newdata$het.scale*newdata$patch
    #browser()
    p_mat[i, ] <- inv_logit(eta)
  }
  #browser()
  # extract the median and credible intervals of predictions
  p_med <- apply(p_mat,2,median)
  p_lower <- apply(p_mat,2,quantile,0.025)
  p_upper <- apply(p_mat,2,quantile,0.975)
  
  # add the medians and CIs to data frame
  newdata$pmed <- p_med
  newdata$plower <- p_lower
  newdata$pupper <- p_upper
  
  # separate by neighbor identity because we want to see relationships for conspecific/heterospecific when other is 0
  only.con <- newdata[newdata$het.raw==0,]
  only.het <- newdata[newdata$con.raw==0,]
  
  # average over patch identity - patch was in stan model but don't care about the relationship in fig 3
  only.con <- only.con %>% group_by(elevation,con.raw) %>% summarize(pmed = mean(pmed),
                                                                     plower=mean(plower),
                                                                     pupper=mean(pupper))
  only.het <- only.het %>% group_by(elevation,het.raw) %>% summarize(pmed = mean(pmed),
                                                                     plower=mean(plower),
                                                                     pupper=mean(pupper))
  #browser()
  # formatting and combining into one data frame
  only.con <- rename(only.con, num.nb = con.raw)
  only.het <- rename(only.het,num.nb=het.raw)
  only.con$neighbors <- "conspecific"
  only.het$neighbors <- "heterospecific"
  ypred <- rbind(only.con,only.het)
  
  return(ypred)
}

predict.srl <- function(stanfit,srl.data){
  
  # extract the posterior samples for all betas
  post <- extract(stanfit)
  beta_elev <- post$beta[,1]
  beta_p <- post$beta[,2]
  beta_con <- post$beta[,3]
  beta_het <- post$beta[,4]
  beta_elxcon <- post$beta[,5]
  beta_elxhet <- post$beta[,6]
  beta_pxcon <- post$beta[,7]
  beta_pxhet <- post$beta[,8]
  
  # find the mean and standard deviation of raw data (I scaled data in stan model)
  mean.raw.con <- mean(srl.data$n.con.nb)
  mean.raw.het <- mean(srl.data$n.het.nb)
  sd.raw.con <- sd(srl.data$n.con.nb)
  sd.raw.het <- sd(srl.data$n.het.nb)
  
  # create new y data frame
  newdata <- expand.grid(
    con.raw = 0:30,
    het.raw = 0:30,
    elevation=c(0,1),
    patch=c(0,1)
  )
  
  # scale the number of neighbors, like I did for stan model
  newdata$con.scale <- (newdata$con.raw - mean.raw.con) / sd.raw.con
  newdata$het.scale <- (newdata$het.raw - mean.raw.het) / sd.raw.het
  
  # set up prediction matrix & prepare for loop
  inv_logit <- function(x){1/(1+exp(-x))}
  ndraws <- length(post$alpha)
  npoints <- nrow(newdata)
  y_mat <- matrix(NA,ndraws,npoints)
  
  # pull beta posterior samples and use new data to estimate over new neighbor abundances
  for(i in 1:ndraws){
    eta <- post$alpha[i] + beta_elev[i]*newdata$elevation + beta_p[i]*newdata$patch + beta_con[i]*newdata$con.scale + beta_het[i]*newdata$het.scale +
      beta_elxcon[i]*newdata$con.scale*newdata$elevation + beta_elxhet[i]*newdata$het.scale*newdata$elevation+
      beta_pxcon[i]*newdata$con.scale*newdata$patch + beta_pxhet[i]*newdata$het.scale*newdata$patch
    #browser()
    y_mat[i, ] <- exp(eta)
  }
  #browser()
  # extract the median and credible intervals of predictions
  y_med <- apply(y_mat,2,median)
  y_lower <- apply(y_mat,2,quantile,0.025)
  y_upper <- apply(y_mat,2,quantile,0.975)
  
  # add the medians and CIs to data frame
  newdata$ymed <- y_med
  newdata$ylower <- y_lower
  newdata$yupper <- y_upper
  
  # separate by neighbor identity because we want to see relationships for conspecific/heterospecific when other is 0
  only.con <- newdata[newdata$het.raw==0,]
  only.het <- newdata[newdata$con.raw==0,]
  
  # average over elevation identity - elevation was in stan model but don't care about the relationship in fig 4
  only.con <- only.con %>% group_by(patch,con.raw) %>% summarize(ymed = mean(ymed),
                                                                     ylower=mean(ylower),
                                                                     yupper=mean(yupper))
  only.het <- only.het %>% group_by(patch,het.raw) %>% summarize(ymed = mean(ymed),
                                                                     ylower=mean(ylower),
                                                                     yupper=mean(yupper))
  #browser()
  # formatting and combining into one data frame
  only.con <- rename(only.con, num.nb = con.raw)
  only.het <- rename(only.het,num.nb=het.raw)
  only.con$neighbors <- "conspecific"
  only.het$neighbors <- "heterospecific"
  ypred <- rbind(only.con,only.het)
  
  return(ypred)
}
