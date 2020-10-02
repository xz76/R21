library(survival)
library(boot)

simulate <- function(n=200, rate=2, c.rate=1, prob.z = 0.4, hr.Z=1.5, 
                     indep_censor = TRUE) {
  Z <- rbinom(n=n, size=1, prob=prob.z)
  rate <- Z*rate*hr.Z + (1-Z)*rate
  if (!indep_censor) {
    c.rate <- Z*c.rate*hr.Z + (1-Z)*c.rate
  }
  T <- rexp(n=n, rate=rate)
  U <- rexp(n=n, rate=c.rate)
  
  X <- mapply(min, T, U)
  # True survival under exponential distribution is exp(-rate*t)
  delta <- 1*(T<=U)
  dat <- as.data.frame(cbind(X, delta, Z))
  return(dat)
}

#IPCW esitmator for independent right censoring
Sipw_t <- function(times=1, X, delta){
  #Estimate censoring distribution 
  c.delta <- 1 - delta
  fit <- survfit(Surv(X, c.delta) ~ 1)
  #Get the value that corresponds to the maximum t that is <= "time"
  element <- function(t){
    max(1,max(1:length(fit$time)*(fit$time<=t)))
  }
  elements<-sapply(times,element)
  G_t <- fit$surv[elements]
  P_n <- function(t){
    mean(X>t)
  }
  
  surv <- sapply(times,P_n)/G_t
  se <- sqrt(surv*(1-surv)/length(X))
  if(length(surv)==1){
    res<-c(surv, se)
  } else {
    res<-rbind(surv, se)
  }
  return(res)
}

#IPCW esitmator for dependent right censoring (through Z)
Sipw_t_z <- function(times=1, X, delta, Z){
  #Estimate censoring distribution 
  c.delta <- 1 - delta
  fit <- coxph(Surv(X, c.delta) ~ Z)
  beta <- coef(fit)
  fit <- basehaz(fit, centered = FALSE)
  
  #Get the value that corresponds to the maximum t that is <= "time"
  element <- function(t){
    max(1,max(1:length(fit$time)*(fit$time<=t)))
  }
 
  P_n <- function(t){
    G_t <- exp(-fit$hazard[element(t)]%*%t(exp(beta*Z)))
    rowMeans((X>t)/G_t)
  }
  
  surv <- sapply(times,P_n)
  se <- sqrt(surv*(1-surv)/length(X))
  if(length(surv)==1){
    res<-c(surv, se)
  } else {
    res<-rbind(surv, se)
  }
  return(res)
}

#Simulation function for hte IPCW estimator for censoring
run.sim2 <- function(N=1000, n=200, rate=2, times=c(0.5, 1, 2), seed=12345, 
                     prob.z=0.4, hr.Z=1.5, indep_censor = TRUE){
  set.seed(seed)
  
  S_n <- matrix(NA, nrow=N, ncol=length(times))
  se_n <- matrix(NA, nrow=N, ncol=length(times))
  ecp <- matrix(NA, nrow=N, ncol=length(times))
  S_0 <- exp(-rate*times)*(1-prob.z) + exp(-rate*hr.Z*times)*prob.z #truth
  for(i in 1:N){
    dat <- simulate(n=n, rate=rate, prob.z = prob.z, hr.Z=hr.Z,
                    indep_censor = indep_censor)
    res <- Sipw_t_z(times=times, X=dat$X, delta=dat$delta, Z=dat$Z)
    S_n[i,] <- res[1,]
    se_n[i,] <- res[2,]
    ecp[i,] <- 1*(S_n[i,] - 1.96*se_n[i,] <= S_0 &
                    S_n[i,] + 1.96*se_n[i,] >= S_0)
  }
  
  
  bias <- colMeans(S_n) - S_0
  MCSD <- apply(S_n, 2, sd)
  ESE <- colMeans(se_n)
  
  CP <- colMeans(ecp)
  
  res <- as.data.frame(cbind(bias, MCSD, ESE, CP))
  colnames(res) <- c("bias", "MCSD", "ESE", "CP")
  rownames(res) <- times
  
  return(res)
}

#################### Bootstrap ######################################
#####################################################################
## Surv Estimation.
Sipw_surv <- function(dat, d = seq_len(nrow(dat)), times=1) {
  X <- dat$X[d]
  delta <- dat$delta[d]
  Z <- dat$Z[d]
  #Estimate censoring distribution 
  c.delta <- 1 - delta
  # fit <- tryCatch(coxph(Surv(X, c.delta) ~ Z), 
  #                 warning = function(w) { w })
  # if (inherits(fit, "warning")) {
  #   browser()
  # }
  fit <- coxph(Surv(X, c.delta) ~ Z)
  beta <- coef(fit)
  fit <- basehaz(fit, centered = FALSE)
  
  #Get the value that corresponds to the maximum t that is <= "time"
  element <- function(t){
    max(1,max(1:length(fit$time)*(fit$time<=t)))
  }
  
  P_n <- function(t){
    G_t <- exp(-fit$hazard[element(t)]%*%t(exp(beta*Z)))
    rowMeans((X>t)/G_t)
  }
  surv <- sapply(times,P_n)
  return(surv)
}

dat <- simulate()
## Bootstrap SE.
boot_se <- function(dat, times = 1, R = 1000) {
  res <- boot(dat, Sipw_surv, R = R, times = times, strata = dat$delta)
  apply(res$t, 2, sd)
}

Sipw_boot <- function(times = 1, X, delta, Z, R = 1000){
  dat <- data.frame(X = X, delta = delta, Z = Z)
  surv <- Sipw_surv(dat, times = times)
  se <- boot_se(dat, times = times, R = R)
  cbind(times, surv, se)
}

run.sim3 <- function(N=1000, n=200, rate=2, times=c(0.5, 1, 2), seed=12345, 
                     prob.z=0.4, hr.Z=1.5, indep_censor = TRUE, R = 1000){
  set.seed(seed)
  
  S_n <- matrix(NA, nrow=N, ncol=length(times))
  se_n <- matrix(NA, nrow=N, ncol=length(times))
  ecp <- matrix(NA, nrow=N, ncol=length(times))
  S_0 <- exp(-rate*times)*(1-prob.z) + exp(-rate*hr.Z*times)*prob.z #truth
  for(i in 1:N){
    dat <- simulate(n=n, rate=rate, prob.z = prob.z, hr.Z=hr.Z,
                    indep_censor = indep_censor)
    res <- Sipw_boot(times=times, X=dat$X, delta=dat$delta, Z=dat$Z, R = R)
    S_n[i,] <- res[, 2]
    se_n[i,] <- res[, 3]
    ecp[i,] <- 1*(S_n[i,] - 1.96*se_n[i,] <= S_0 &
                    S_n[i,] + 1.96*se_n[i,] >= S_0)
  }
  bias <- colMeans(S_n) - S_0
  MCSD <- apply(S_n, 2, sd)
  ESE <- colMeans(se_n)
  
  CP <- colMeans(ecp)
  
  res <- as.data.frame(cbind(bias, MCSD, ESE, CP))
  colnames(res) <- c("bias", "MCSD", "ESE", "CP")
  rownames(res) <- times
  
  return(res)
}

#Biased
#run.sim1(N=1000, n=200, rate=1, times = c(0.5, 1, 1.5))

#Censoring adjusted estimator (but SEs are invalid since the y ingore the variability in the KM esitmator of censoring)
org_result <- run.sim2(N=100, n=200, rate=1, times = c(0.5, 1, 1.5), 
                       indep_censor = TRUE)

boot_result2 <- run.sim3(N=100, n=200, rate=1, times = c(0.5, 1, 1.5), 
                        indep_censor = TRUE, R = 1000)
saveRDS(org_result, file = "original_result.rds")
saveRDS(boot_result2, file = "boot_result2.rds")

