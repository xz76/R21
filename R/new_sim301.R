# Estimator of state occupation probabilities with:
# - Dependent right censoring: X(t) <- Z -> CC

library(truncdist)
library(survival)
library(Matrix)
library(foreach)
library(doMC)

doMC::registerDoMC(cores = parallel::detectCores())

# Simulate right-censored processes from the illness-death model with recovery
simulate <- function(n=1000, a12=0.5, a13=1, a21=0.75 ,a23=0.5,
                     p12=1, p13=1, p21=1, p23=1, c.rate=1, b.Z.cens=0.5,
                     p.Z=0.4, b.Z=0.5, gamma_0=-1, gamma_1=1, scale0 = 1) {
  # In this function when p12=p13=p21=p23=1 we have a homogeneous process
  # (i.e. time-constant intensities)

  outdata <- NULL
  for(i in 1:n){
    #Simulate covariate that is associated with X(t), right censoring, and missingness
    Z <- rbinom(n=1, size=1, prob=p.Z)

    #Probability of missingness
    p.miss <- exp(gamma_0 + gamma_1*Z)/(1+exp(gamma_0 + gamma_1*Z))
    #R <- rbinom(n=1, size=1, prob=p.miss)
    R=1
    if(R==1){

      CC <- rexp(1, rate=c.rate*exp(b.Z.cens*Z))

      t12 <- rweibull(1, shape=p12, scale=(scale0/(a12*exp(b.Z*Z))))
      t13 <- rweibull(1, shape=p13, scale=(scale0/(a13*exp(b.Z*Z))))
      stop <- 1*((t13 < t12) | (CC < t12)) #Stop if death or censoring first
      x1 <- min(c(t12, t13)) # Exit time from state 1 (ignoring censoring)
      if (stop==0){
        t1 <- 0
        t2 <- x1
        s1 <- 1
        s2 <- 2
        while (stop==0){
          if(s2[length(s2)]==2){
            t23 <- rtrunc(1, spec="weibull", a=t2[length(t2)], b=Inf,
                          shape=p23, scale=(scale0/(a23*exp(b.Z*Z))))
            t21 <- rtrunc(1, spec="weibull", a=t2[length(t2)], b=Inf,
                          shape=p21, scale=(scale0/(a21*exp(b.Z*Z))))

            stop <- 1*((t23 < t21) | (CC < t21)) #Stop if death or censoring first
            x2 <- min(c(t21, t23)) # Exit time from state 2 (ignoring censoring)

            t1 <- c(t1, x1)
            t2 <- c(t2, min(x2, CC))
            s1 <- c(s1, 2)
            s2 <- c(s2, ifelse(x2 <= CC, ifelse(t21 <= t23, 1, 3), 2))
          } else {
            t12 <- rtrunc(1, spec="weibull", a=t2[length(t2)], b=Inf,
                          shape=p12, scale=(scale0/(a12*exp(b.Z*Z))))
            t13 <- rtrunc(1, spec="weibull", a=t2[length(t2)], b=Inf,
                          shape=p13, scale=(scale0/(a13*exp(b.Z*Z))))

            stop <- 1*((t13 < t12) | (CC < t12)) #Stop if death or censoring first
            x1 <- min(c(t12, t13)) # Exit time from state 1 (ignoring censoring)

            t1 <- c(t1, x2)
            t2 <- c(t2, min(x1, CC))
            s1 <- c(s1, 1)
            s2 <- c(s2, ifelse(x1 <= CC, ifelse(t12 <= t13, 2, 3), 1))
          }
        }
      } else {
        t1 <- 0
        t2 <- min(x1, CC)
        s1 <- 1
        s2 <- ifelse(x1 < CC, 3, 1)
      }
    } else {
      t1 <- 0
      t2 <- NA
      s1 <- 0
      s2 <- NA
    }

    id <- rep(i, length(t1))
    Z <- rep(Z, length(t1))
    R <- rep(R, length(t1))
    outdata <- rbind(outdata, data.frame(id, t1, t2, s1, s2, Z, R))
  }
  # t1: start time of the interval
  # t2: end time of the interval
  # s1: state at the start time of the interval
  # s2: state at the end time of the interval
  # if s1==s2 then right censoring occured at t2 at this interval
  outdata <- outdata[order(outdata$id, outdata$t1, outdata$t2),]
  # outdata$t1 <- 1000 * outdata$t1
  # outdata$t2 <- 1000 * outdata$t2
  return(outdata)
}


# transmatrix function
trans <- function(nstate, state_names, from, to) {
  if (missing(nstate) && missing(state_names))
    stop("One of 'nstate' and 'state_names' has to be specified.")
  if (missing(state_names)) {
    state_names <- as.character(seq_len(nstate))
  } else {
    state_names <- unique(state_names)
    nstate <- length(state_names)
  }
  if (length(from) != length(to))
    stop("The length of 'from' and 'to' must be the same.")
  if (is.character(from)) {
    from <- match(from, state_names)
  } else {
    from <- as.integer(from)
  }
  if (is.character(to)) {
    to <- match(to, state_names)
  } else {
    to <- as.integer(to)
  }
  mat <- matrix(FALSE, ncol = nstate, nrow = nstate)
  dimnames(mat) <- list(state_names, state_names)
  mat[cbind(from, to)] <- TRUE
  mat
}


### New sop function
sop_tnew <- function(data, tau=NULL, S, T_c, ipw=0, trans,
                     times= c(0.1, 0.2)){
    if(is.null(tau)){
      tau <- max(data$t2)
    }
    CTI <- list()
    counter <- 1
    for(h in T_c){
      for(j in S[trans[h,]]){
        data_h <- data[data$s1==h,]
        data_h$delta <- 1*(data_h$s2==j)
        if(ipw==0){
          fit <- coxph(Surv(t1,t2,delta,type="counting")~1, data=data_h,
                       control = coxph.control(timefix = FALSE))
        } else {
          fit <- coxph(Surv(t1,t2,delta,type="counting")~1,
                       weight=weightVL, data=data_h,
                       control = coxph.control(timefix = FALSE))
        }
        A <- basehaz(fit, centered=FALSE)
        A_t<-stepfun(A$time,c(0,A$hazard))
        CTI[[counter]] <- A_t
        if(counter==1){
          pointer <- c(h,j,counter)
        } else {
          pointer <- rbind(pointer,c(h,j,counter))
        }
        counter <- counter + 1
      }
    }

    tt<-sort(unique(data[data$s2 != data$s1,"t2"]))
    tt <- tt[tt<=tau]

    dA <- sapply(seq_along(CTI), function(i) {
      diff(c(0, CTI[[i]](tt)), lag = 1)
    })

    ttrans <- t(trans)
    mat0 <- matrix(0, nrow = nrow(trans), ncol = ncol(trans))
    mat_list <- lapply(seq_len(nrow(dA)), function(i) {
      out <- mat0
      out[which(ttrans)] <- dA[i, ]
      out <- t(out)
      ## compute diagonal elements
      diag(out) <- - rowSums(out)
      out + diag(nrow(trans))
    })

    #Calculate the product integral estimator
    P_n<-Reduce("%*%",  mat_list, accumulate = TRUE)
    p_0 <- sapply(S, function(i){sum(data[data$t1==0,"s1"] == i)/nrow(data[data$t1==0,])})
    p_n <- matrix(NA,nrow=length(tt),ncol= nrow(trans))
    for(i in 1:length(tt)){
      p_n[i,] <- p_0%*%P_n[[i]]
    }
    p_n <- rbind(c(0,p_0),
                 cbind(tt,p_n))
    p_n <- as.data.frame(p_n)
    colnames(p_n) <- c("t",paste("p",S,sep=""))
    rownames(p_n) <- 1:(length(tt)+1)
    if (!is.null(times)){
      p_nt <- sapply(times, function(t){tail(p_n[which(p_n$t < t), ], 1)})
      colnames(p_nt) <- times
      res <- t(p_nt)
      out <- do.call(c, res)
      attributes(out) <- attributes(res)
    }else{
      out <- p_n
    }
   
    return(out)
}

se_bootnew <- function(data, times=1, j, nboot=100,
                    tau = NULL, S = 1:3, T_c = 1:2, ipw = 0, trans){
  id <- sort(unique(data$id))
  theta <- NULL
  for(b in 1:nboot){
    bdat <- NULL
    bid <- sample(id, replace=TRUE)
    for(i in 1:length(bid)){
      bdat <- rbind(bdat, data[data$id==bid[i],])
    }
    if (all(bdat$s1 == bdat$s2)){
      browser()
      next
    }
      res <- sop_tnew(data = bdat, tau=tau, S = S, T_c = T_c, ipw= ipw,
                      trans = trans, times= times)
      theta <- rbind(theta, res)
  }
  dat <- lapply(S, function(i){
    idx <- seq(i, nrow(theta), by = 3)
    theta[idx,]
  })
  b.se <- lapply(dat, function(data){
   apply(data, 2, sd)
  })
  b.se <- do.call(rbind, b.se)
  b.se <- b.se[,-1]
  return(b.se)
}


# Monte-Carlo calculation of the true P_j(t)
# n should be large in order to have an accurate approximation of
# P_j(t) = E[I{X(t)=j}]
true.P_j <- function(n=10000, a12=0.5, a13=1, a21=0.25 ,a23=0.5,
                     p12=1, p13=1, p21=1, p23=1,
                     p.Z=0.4, b.Z=0.5, j=2, times, scale0 = 0.01) {
    # In this function when p12=p13=p21=p23=1 we have a homogeneous process
    # (i.e. time-constant intensities)
    if(p12==1 & p13==1 & p21==1 & p23==1){
      a_0 <- rbind(c(-(a12+a13), a12, a13),
                   c(a21, -(a21+a23), a23),
                   c(0, 0, 0))
      P_t <- function(t){
        expm(a_0*t)[1, j]
      }
      P0 <- sapply(times, P_t)
      P_t <- function(t){
        expm(a_0*exp(b.Z)*t)[1, j]
      }
      P1 <- sapply(times, P_t)

      #True (marginal on Z) state occupation probability
      P_0 <- p.Z*P1 + (1 - p.Z)*P0

    }else{
      outdata <- NULL
      for(i in 1:n){
        #Simulate covariate that is associated with X(t), right censoring, and missingness
        Z <- rbinom(n=1, size=1, prob=p.Z)
        t12 <- rweibull(1, shape=p12, scale=(scale0/(a12*exp(b.Z*Z))))
        t13 <- rweibull(1, shape=p13, scale=(scale0/(a13*exp(b.Z*Z))))
        stop <- 1*(t13 < t12) #Stop if death first
        x1 <- min(c(t12, t13)) # Exit time from state 1
        if (stop==0){
          t1 <- 0
          t2 <- x1
          s1 <- 1
          s2 <- 2
          while (stop==0){
            if(s2[length(s2)]==2){
              t23 <- rtrunc(1, spec="weibull", a=t2[length(t2)], b=Inf,
                            shape=p23, scale=(scale0/(a23*exp(b.Z*Z))))
              t21 <- rtrunc(1, spec="weibull", a=t2[length(t2)], b=Inf,
                            shape=p21, scale=(scale0/(a21*exp(b.Z*Z))))

              stop <- 1*(t23 < t21) #Stop if death first
              x2 <- min(c(t21, t23)) # Exit time from state 2

              t1 <- c(t1, x1)
              t2 <- c(t2, x2)
              s1 <- c(s1, 2)
              s2 <- c(s2, ifelse(t21 <= t23, 1, 3))
            } else {
              t12 <- rtrunc(1, spec="weibull", a=t2[length(t2)], b=Inf,
                            shape=p12, scale=(scale0/(a12*exp(b.Z*Z))))
              t13 <- rtrunc(1, spec="weibull", a=t2[length(t2)], b=Inf,
                            shape=p13, scale=(scale0/(a13*exp(b.Z*Z))))

              stop <- 1*(t13 < t12) #Stop if death first
              x1 <- min(c(t12, t13)) # Exit time from state 1

              t1 <- c(t1, x2)
              t2 <- c(t2, x1)
              s1 <- c(s1, 1)
              s2 <- c(s2, ifelse(t12 <= t13, 2, 3))
            }
          }
        } else {
          t1 <- 0
          t2 <- x1
          s1 <- 1
          s2 <- 3
        }
        id <- rep(i, length(t1))
        Z <- rep(Z, length(t1))
        outdata <- rbind(outdata, data.frame(id, t1, t2, s1, s2, Z))
      }
      # t1: start time of the interval
      # t2: end time of the interval
      # s1: state at the start time of the interval
      # s2: state at the end time of the interval
      # if s1==s2 then right censoring occurred at t2 at this interval
      outdata <- outdata[order(outdata$id, outdata$t1, outdata$t2),]

      #Average state occupation probability
      P_t <- function(t){
        Ind_j <- (j<3)*(outdata$s1==j & outdata$t1<=t & outdata$t2>t) +
          (j==3)*(outdata$s2==j & outdata$t2<=t)
        sum(Ind_j)/n
      }
      P_0 <- sapply(times, P_t)
    }
  return(P_0)
}



#Simulation function for the IPCW estimator for dependent censoring
run.sim <- function(N=1000, n=1000, times = c(0.5, 1, 1.5), state=2,
                     a12=0.5, a13=1, a21=0.75 ,a23=0.5, seed = 123,
                     p12=1, p13=1, p21=1, p23=1, c.rate=1, b.Z.cens=0.5,
                     p.Z=0.4, b.Z=0.5, gamma_0=-1, gamma_1=1, nboot=100,
                    trace=TRUE,  tau=NULL, S, T_c , ipw, trans, P_0 ){
  set.seed(seed)
  out <- foreach(i = seq_len(N)) %dopar%
    ({ 
    dat <- simulate(n=n, a12=a12, a13=a13, a21=a21, a23=a23,
                         p12=p12, p13=p13, p21=p21, p23=p23, c.rate=c.rate,
                         b.Z.cens=b.Z.cens, p.Z=p.Z, b.Z=b.Z,
                         gamma_0=gamma_0, gamma_1=gamma_1)
    res <- sop_tnew(dat, tau=tau, S = S, T_c = T_c, ipw= ipw,
                    trans = trans, times= times)
    se <- se_bootnew(dat, times = times, nboot = nboot, tau = tau, S = S,
                     T_c = T_c, ipw = ipw, trans = trans, j = state)
    se_n <- se[, state]
    P_n <- res[, state +1]
    ecp <- 1*(P_n - 1.96*se_n <= P_0 &
              P_n + 1.96*se_n >= P_0)

    list(P_n = P_n, se_n = se_n, ecp = ecp)
    })
  P_n <- do.call(rbind, lapply(out, function(a) a$P_n))
  se_n <- do.call(rbind, lapply(out, function(a) a$se_n))
  ecp <- do.call(rbind, lapply(out, function(a) a$ecp))

  bias <- 100*(colMeans(P_n) - P_0)/P_0
  MCSD <- apply(P_n, 2, sd)
  ESE <- colMeans(se_n)
  CP <- colMeans(ecp)

  res <- as.data.frame(cbind(bias, MCSD, ESE, CP))
  colnames(res) <- c("bias", "MCSD", "ESE", "CP")
  rownames(res) <- times
  return(res)
}
#Censoring adjusted estimator
set.seed(12345)
tmatrix <- trans(state_names = c("health", "illness", "death"),from = c(1, 1, 1, 2, 2),
                     to = c(2, 2, 3, 3, 1))
# P_0 <- true.P_j( n=10000, times = c(0.5, 1, 1.5), j=2,
#                  a12=0.5, a13=1, a21=0.75 ,a23=0.5,
#                  p12=1, p13=1, p21=1, p23=1)

# res_200 <- run.sim(N=1000, n=200, times = c(0.5, 1, 1.5), trace=FALSE,  tau=NULL, S = 1:3,
#         T_c = 1:2, ipw=0, trans = tmatrix, P_0 = P_0)
# 
# res_400 <- run.sim(N=1000, n=400, times = c(0.5, 1, 1.5), trace=FALSE,  tau=NULL, S = 1:3,
#         T_c = 1:2, ipw=0, trans = tmatrix, P_0 = P_0)
# res_800 <- run.sim(N=1000, n=800, times = c(0.5, 1, 1.5), trace=FALSE,  tau=NULL, S = 1:3,
#                    T_c = 1:2, ipw=0, trans = tmatrix, P_0 = P_0)
#saveRDS(res_200, file = "res_200.rds")
# saveRDS(res_400, file = "res_400.rds")
# saveRDS(res_800, file = "res_800.rds")

## Event Numbers
# event_rate <- function(N, n, times = c(0.5, 1, 1.5), j=2,
#                   a12=0.5, a13=1, a21=0.75 ,a23=0.5,
#                   p12=0.75, p13=1, p21=1, p23=1, from = 1, to = 2){
#   
#   res <- replicate(N, { sapply(times, function(t){
#     tmp = simulate(n = n, p12 = p12)
#     sum(tmp$s1 == from & tmp$s2 == to & tmp$t2 <= t)
#   })
#   })
#   return(apply(res, 1, mean))
# }
# e.num200 <- event_rate(N = 1000, n = 200, p12 = 0.75)
# e.num400 <- event_rate(N = 1000, n = 400, p12 = 0.75)
# e.num800 <- event_rate(N = 1000, n = 800, p12 = 0.75)

P_1 <- true.P_j( n=10000, times = c(0.5, 1, 1.5), j=2,
                 a12=0.5, a13=1, a21=0.75 ,a23=0.5,
                 p12=0.55, p13=1, p21=1, p23=1, scale0 = 1)


# try200 <- run.sim(N=10, n=200, times = c(0.5, 1, 1.5), trace=FALSE,  tau=NULL, S = 1:3,
#                   T_c = 1:2, ipw=0, trans = tmatrix, P_0 = P_1, p12=0.55)

nonhomo_200 <- run.sim(N=500, n=200, times = c(0.5, 1, 1.5), trace=FALSE,  tau=NULL, S = 1:3,
                       T_c = 1:2, ipw=0, trans = tmatrix, P_0 = P_1,  p12=0.55)
saveRDS(nonhomo_200, file = "nonhomo_200.rds")

nonhomo_400 <- run.sim(N=500, n=400, times = c(0.5, 1, 1.5), trace=FALSE,  tau=NULL, S = 1:3,
                       T_c = 1:2, ipw=0, trans = tmatrix, P_0 = P_1, p12 = 0.55)
saveRDS(nonhomo_400, file = "nonhomo_400.rds")

nonhomo_800 <- run.sim(N=400, n=800, times = c(0.5, 1, 1.5), trace=FALSE,  tau=NULL, S = 1:3,
                       T_c = 1:2, ipw=0, trans = tmatrix, P_0 = P_1, p12=0.55)

saveRDS(nonhomo_800, file = "nonhomo_800.rds")


###########
##### Statistics
# d1 <- simulate(n = 200,  p12=0.5, p23 = 0.7)
# d2 <- simulate(n = 200, p12 = 0.35)
# stats <- function(data, timepoint = seq(0, 10, 1)){
#   state_num <- function(time, state, dat) {
#     nrow(subset(dat, s1 == state & t1 <= time & time < t2))
#   }
#   state1 <- sapply(timepoint, state_num, state = 1, dat = data)
#   state2 <- sapply(timepoint, state_num, state = 2, dat = data)
#   res <- data.frame(timepoint, state1, state2)
#   res
# }
# 
# res$proportion <- percent(res$state1/res$state2)
