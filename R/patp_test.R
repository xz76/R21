## Function to calculate the p-value for the two-sample comparison
## of transition probabilities based on a Kolmogorov-Smirnov-type
## test

simulate <- function(n=1000, a12=0.5, a13=1, a21=0.75 ,a23=0.5,
                     p12=1, p13=1, p21=1, p23=1, c.rate=1, b.Z.cens=0.5,
                     p.Z=0.4, b.Z=0.5, gamma_0=-1, gamma_1=1) {
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

      t12 <- rweibull(1, shape=p12, scale=(1/(a12*exp(b.Z*Z))))
      t13 <- rweibull(1, shape=p13, scale=(1/(a13*exp(b.Z*Z))))
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
                          shape=p23, scale=(1/(a23*exp(b.Z*Z))))
            t21 <- rtrunc(1, spec="weibull", a=t2[length(t2)], b=Inf,
                          shape=p21, scale=(1/(a21*exp(b.Z*Z))))

            stop <- 1*((t23 < t21) | (CC < t21)) #Stop if death or censoring first
            x2 <- min(c(t21, t23)) # Exit time from state 2 (ignoring censoring)

            t1 <- c(t1, x1)
            t2 <- c(t2, min(x2, CC))
            s1 <- c(s1, 2)
            s2 <- c(s2, ifelse(x2 <= CC, ifelse(t21 <= t23, 1, 3), 2))
          } else {
            t12 <- rtrunc(1, spec="weibull", a=t2[length(t2)], b=Inf,
                          shape=p12, scale=(1/(a12*exp(b.Z*Z))))
            t13 <- rtrunc(1, spec="weibull", a=t2[length(t2)], b=Inf,
                          shape=p13, scale=(1/(a13*exp(b.Z*Z))))

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

#######################################
#### Add group and cid information#####
#######################################
sim_group <- function(data, cid_num, trans){
  l <- length(unique(data$id))
  g <- rbinom(l, 1, 0.5)
  q <- quantile(l)
  rand_cid <- function(cid_num, l, sd = 1) {
    vec <- rnorm(cid_num - 1, l/(cid_num), sd)
    vec <-  floor(vec)
    vec[cid_num] <- l - sum(vec)
    vec
  }
  cid_seq <- rand_cid(cid_num = cid_num, l = l ,sd = 1)
  ## Add group information
  res <- do.call(rbind, by(data, data$id, function(sub_dat){
    sub_dat$group <- rep(g[sub_dat$id[1]], nrow(sub_dat))
    return(sub_dat)
  }))
  cid <- NULL
  ## Add cid information
  for(i in 1 : length(cid_seq)){
    tmp <- replicate(cid_seq[i], i)
    cid <- c(cid,tmp)
  }
  out <- do.call(rbind, by(res, res$id, function(sub_dat){
    sub_dat$cid <- rep(cid[sub_dat$id[1]], nrow(sub_dat))
    return(sub_dat)
  }))
  rownames(out) <- NULL
  return(out)
}

tmat1 <- transMat(x = list(c(2, 3), c(3), c()), names = c("Health", "Illness", "Death"))
tmp <- simulate(100)
tmp <- sim_group(tmp, cid = 3)
data <- tmp

data$trans <- with(data, tmat1[cbind(s1, s2)])

patp_test <- function(data, tmat = tmat1, cid, id, group, h, j, s=0,
                      weighted=FALSE, LMAJ=FALSE, B=1000){
  # microbenchmark::microbenchmark(aggregate(data[,cid], by=list(data[,id]),
  #                                          FUN=sd, na.rm=TRUE)$x,
  #                                tapply(data[, cid], data[, id], sd, na.rm = TRUE))
  # check.ic <- aggregate(data[,cid], by=list(data[,id]),
  #                       FUN=sd, na.rm=TRUE)$x
  ## tapply faster than aggregate function.

  check.ic <- tapply(data[, cid], data[, id], sd, na.rm = TRUE)
  check.ic <- check.ic[!is.na(check.ic)]
  if(length(check.ic)>0){
    if(max(check.ic)>0){
      stop("Same unit(s) in more than 1 cluster (violation of the independent clusters assumption)")
    }
  }

  if(B<=0){
    stop("Tests cannot be performed based on <=0 bootstrap samples")
  } else if (B<1000){
    warning("It is recommended to use at least 1000 bootstrap samples for two-sample testing")
  }
  n <- length(unique(data[,cid]))
  groups <- unique(data[,group])
  groups <- sort(groups[!is.na(groups)])
  if(length(groups)!=2){
    stop("Number of groups != 2")
  }

  if(weighted==TRUE){
    M0 <- aggregate(rep(1,times=nrow(data)),
                    by = list(data[,cid],
                              data[,group],
                              data[,id]),
                    FUN = mean)
    # with(data, tapply(rep(1,times=nrow(data)),
    #                         paste(cid, group, id, sep = ":"), FUN = mean))
    M <- aggregate(M0$x,
                   by = list(M0$Group.1, M0$Group.2),
                   FUN = sum)
    colnames(M) <- c(cid, group, "clust.size")
    data <- merge(data, M, by=c(cid,group))
    data <- data[order(data[,cid],data[,id]),]
    class(data) <- c("msdata", "data.frame")
  }

  Pt <- list()
  tms <- list()
  for(g in groups){
    if(LMAJ==FALSE){
      c0 <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans),
                  data=data[data[,group]==g,], method = "breslow")

      A0 <- msfit(object = c0, trans = tmat, variance=FALSE)
      if(weighted==TRUE){
        c0 <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans),
                    weights=(1/clust.size), data=data[data[,group]==g,],
                    method = "breslow")
        A.wt <- basehaz(c0, centered=FALSE)
        A.wt$strata <- as.numeric(A.wt$strata)
        for(trn in sort(unique(A.wt$strata))){
          fun <- stepfun(A.wt[A.wt$strata==trn,"time"],
                         c(0,A.wt[A.wt$strata==trn,"hazard"]))
          A0$Haz[A0$Haz$trans==trn,"Haz"] <- fun(A0$Haz[A0$Haz$trans==trn,"time"])
        }
      }
      #######
      P0 <- sop_tnew()
      ######

      P0 <- probtrans(A0, predt = s,
                      variance=FALSE)[[h]][,c("time", paste("pstate", j, sep=""))]
    } else {
      P0 <- LMAJ2(msdata=data[data[,group]==g,], tmat=tmat,
                  id=id, s=s, h=h, j=j, weighted=weighted)
    }
    if(length(Pt)==0){
      Pt[[1]] <- stepfun(P0$time, c(P0[1,paste("pstate", j, sep="")],
                                    P0[,paste("pstate", j, sep="")]))
      tms[[1]] <- P0$time[diff(c(P0[1,paste("pstate", j, sep="")],
                                 P0[,paste("pstate", j, sep="")]), lag=1)!=0]
    } else {
      Pt[[2]] <- stepfun(P0$time, c(P0[1,paste("pstate", j, sep="")],
                                    P0[,paste("pstate", j, sep="")]))
      tms[[2]] <- P0$time[diff(c(P0[1,paste("pstate", j, sep="")],
                                 P0[,paste("pstate", j, sep="")]), lag=1)!=0]
    }
  }

  tms <- sort(unique(c(tms[[1]], tms[[2]])))

  if(nrow(data[data$from==j,])>0){
    tS <- sort(unique(c(data[data$to==j,"from"],j)))
  } else {
    tS <- sort(unique(data[data$to==j,"from"]))
  }

  EY <- NULL
  EY.t <- function(t,dt){
    sum(dt$Tstart<t & dt$Tstop>=t)/n
  }
  for(i in tS){
    for(g in groups){
      dat <- unique(data[data$from==i & data[,group]==g,
                         c(id, "from", "Tstart", "Tstop")])
      if(is.null(EY)){
        EY <- sapply(tms, EY.t, dt=dat)
      } else {
        EY <- cbind(EY, sapply(tms, EY.t, dt=dat))
      }
    }
  }
  Wt <- rowprods(EY)/rowSums(EY)
  tms <- tms[!is.na(Wt)]
  if(length(tms)==0 | max(Wt, na.rm=TRUE)==0){
    stop("Weights NA or 0 for all timepoints")
  }
  Wt <- Wt[!is.na(Wt)]
  D_hat=(Pt[[1]](tms) - Pt[[2]](tms))

  Diff_boot <- msboot(patp_test_b, data=data,
                      id=cid, B=B, verbose=0,
                      tmat=tmat, id2=id, group=group, h=h, j=j, s=s,
                      times=tms, D_hat=D_hat, Wt=Wt,
                      n=n, weighted=weighted, LMAJ=LMAJ)

  KS <- max(abs(sqrt(n)*Wt*D_hat))
  KS.b <- apply(abs(Diff_boot),2,max)
  pval <- mean(KS.b >= KS)
  names(pval) <- "p-value"
  return(pval)
}
