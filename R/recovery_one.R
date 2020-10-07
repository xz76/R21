## read in arguments from bash
in_args <- commandArgs(trailingOnly = TRUE)
data_id <- as.integer(in_args[1])

library(truncdist)
library(mstate)
library(survival)
source("LMAJ2.R")
source("patp_b.R")
source("patp_test_b.R")
source("patp_test.R")
source("patp.R")
# Simulate right-censored processes from the illness-death model with recovery
simulate <- function(n=1000, a12=0.5, a13=1, a21=0.75 ,a23=0.5,
                     p12=1, p13=1, p21=1, p23=1, a=0, b=3) {
  # In this function when p12=p13=p21=p23=1 we have a homogeneous process
  # (i.e. time-constant intensities)
  outdata <- NULL
  for(i in 1:n){
    CC <- runif(1, min=a, max=b)
    t12 <- rweibull(1, shape=p12, scale=(1/a12))
    t13 <- rweibull(1, shape=p13, scale=(1/a13))
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
                        shape=p23, scale=(1/a23))
          t21 <- rtrunc(1, spec="weibull", a=t2[length(t2)], b=Inf,
                        shape=p21, scale=(1/a21))

          stop <- 1*((t23 < t21) | (CC < t21)) #Stop if death or censoring first
          x2 <- min(c(t21, t23)) # Exit time from state 2 (ignoring censoring)

          t1 <- c(t1, x1)
          t2 <- c(t2, min(x2, CC))
          s1 <- c(s1, 2)
          s2 <- c(s2, ifelse(x2 <= CC, ifelse(t21 <= t23, 1, 3), 2))
        } else {
          t12 <- rtrunc(1, spec="weibull", a=t2[length(t2)], b=Inf,
                        shape=p12, scale=(1/a12))
          t13 <- rtrunc(1, spec="weibull", a=t2[length(t2)], b=Inf,
                        shape=p13, scale=(1/a13))

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
    id <- rep(i, length(t1))
    outdata <- rbind(outdata, data.frame(id, t1, t2, s1, s2))
  }
  # t1: start time of the interval
  # t2: end time of the interval
  # s1: state at the start time of the interval
  # s2: state at the end time of the interval
  # if s1==s2 then right censoring occured at t2 at this interval
  outdata <- outdata[order(outdata$id, outdata$t1, outdata$t2),]
  return(outdata)
}

##### Long Form ####
new_long <- function(sub_dat, tmat) {
  mat <- lapply(seq_len(ncol(tmat)), function(i) {
    res <- which(!is.na(tmat[i, ]))
    if (length(res)) {
      out <- cbind(from = i, to = res)
      rownames(out) <- NULL
      out
    } else {
      NULL
    }
  })
  mat <- do.call(rbind, mat)
  from <- mat[, 1]
  tmp_dat <- data.frame(cbind(from = mat[, 1], to = mat[, 2],
                              trans = tmat[mat]))
  
  id <- unique(sub_dat$id)
  is_censor <- sub_dat$s1 == sub_dat$s2

  out <- lapply(seq_along(is_censor), function(i) {
      this_censor <- is_censor[i]
      if (this_censor) {
          res <- data.frame(id = id,
                            Tstart = sub_dat$t1[i],
                            Tstop = sub_dat$t2[i],
                            time = sub_dat$t2[i] - sub_dat$t1[i],
                            from = sub_dat$s1[i])
          res <- merge(res, tmp_dat, by = "from")
          res <- res[, c("id", "from", "to", "Tstart", "Tstop", "time", 
                          "trans")]
          res$status <- 0
      } else {
          res <- data.frame(id = id,
                            Tstart = sub_dat$t1[i],
                            Tstop = sub_dat$t2[i],
                            from = sub_dat$s1[i],
                            to = sub_dat$s2[i],
                            time = sub_dat$t2[i] - sub_dat$t1[i],
                            status = 1)
          res <- merge(res, tmp_dat, by = c("from", "to"))
          res <- res[, c("id", "from", "to", "trans", "Tstart", 
                         "Tstop", "time", "status")]
      }
      res
  })
  do.call(rbind, out)
}

reshape_long <- function(dat) {
  res <- by(data = dat, dat$id, new_long, tmat = tmat)
  res <- do.call(rbind, res)
  class(res) <- c("msdata", "data.frame")
  return(res)
}

true.P_j <- function(n=100000, a12=0.5, a13=1, a21=0.25 ,a23=0.5,
                     p12=1, p13=1, p21=1, p23=1, 
                     j=2, times) {
  # In this function when p12=p13=p21=p23=1 we have a homogeneous process
  # (i.e. time-constant intensities)
  outdata <- NULL
  for(i in 1:n){
    #There is no right-censoring here so that the sample average
    #estimates consistently the expectation E[I{X(t)=j}]
    t12 <- rweibull(1, shape=p12, scale=(1/a12))
    t13 <- rweibull(1, shape=p13, scale=(1/a13))
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
                        shape=p23, scale=(1/a23))
          t21 <- rtrunc(1, spec="weibull", a=t2[length(t2)], b=Inf,
                        shape=p21, scale=(1/a21))
          
          stop <- 1*(t23 < t21) #Stop if death first
          x2 <- min(c(t21, t23)) # Exit time from state 2 
          
          t1 <- c(t1, x1)
          t2 <- c(t2, x2)
          s1 <- c(s1, 2)
          s2 <- c(s2, ifelse(t21 <= t23, 1, 3))
        } else {
          t12 <- rtrunc(1, spec="weibull", a=t2[length(t2)], b=Inf,
                        shape=p12, scale=(1/a12))
          t13 <- rtrunc(1, spec="weibull", a=t2[length(t2)], b=Inf,
                        shape=p13, scale=(1/a13))
          
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
    outdata <- rbind(outdata, data.frame(id, t1, t2, s1, s2))
  }
  # t1: start time of the interval
  # t2: end time of the interval
  # s1: state at the start time of the interval
  # s2: state at the end time of the interval
  # if s1==s2 then right censoring occured at t2 at this interval
  outdata <- outdata[order(outdata$id, outdata$t1, outdata$t2),]
  
  P_t <- function(t){
    Ind_j <- 1*(outdata$s1==j & outdata$t1<=t & outdata$t2>t)
    sum(Ind_j)/n
  }
  P_0 <- sapply(times, P_t)
  return(P_0)
}

# true.P_j(n=1000, j=2, times=c(0,0.01,0.1,.5, 1, 2))

###########
###########

run.sim <- function(dat_long, tmat,
                    n=10000, a12=0.5,a13=1, a21=0.25 ,a23=0.5,
                    p12=1, p13=1, p21=1, p23=1, 
                    j=2, times = c(0.5, 1, 2))
{
  P_0 <-true.P_j(n, a12, a13, a21, a23, p12, p13, p21, p23, j, times)

  P <- patp(data = dat_long, tmat = tmat, cid = "id", id = "id",
              h = 1, j = 2, s = 0, B = 1000)
  
  res <- sapply(times, function(t){P[which.min(abs(P$time - t)),]})
  P_n = unlist(res[2, ])
  se_n = unlist(res[3,])
  ecp =  as.integer(P_n - 1.96 * se_n <= P_0 &
                      P_n + 1.96 * se_n >= P_0)
  res <- data.frame(P_n, se_n, ecp )
  rownames(res) <- times
  return(res)
}
set.seed(123 + data_id)

tmat <- transMat(x = list(c(2, 3), c(1, 3), c()),    
                 names = c("Health", "Illness", "Death"))

temp <- simulate(1000)
res <- reshape_long(temp)
out <- run.sim(res, tmat,
        n=10000, a12=0.5,a13=1, a21=0.25 ,a23=0.5,
        p12=1, p13=1, p21=1, p23=1, 
        j=2, times = c(0.5, 1, 2))

## save the results
out_dir <- "recovery"
if (! dir.exists(out_dir)) {
  dir.create(out_dir)
}
saveRDS(out, file = file.path(out_dir, sprintf("one_sim_%d.rds", data_id)))
