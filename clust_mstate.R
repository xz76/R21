
library(mstate)
library(Rfast)


LMAJ2 <- function (msdata, tmat, id, s, h, j, weighted){
  if (is.null(tmat)) 
    stop("msdata object should have a \"trans\" attribute")
  K <- nrow(tmat)
  if (any(is.na(match(h, 1:K)))) 
    stop("h should be subset of 1:K with K number of states")
  attr(msdata, "trans") <- tmat
  xss <- xsect(msdata, s)
  infrom <- xss[xss$state %in% h,id]
  msdatas <- cutLMms(msdata, LM = s)
  msdatasfrom <- msdatas[msdatas[,id] %in% infrom, ]
  c0 <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), 
              data = msdatasfrom)
  A0 <- msfit(c0, trans = tmat, variance=FALSE)
  if(weighted==TRUE){
    c0 <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), 
                weights=(1/clust.size), data = msdatasfrom)
    A.wt <- basehaz(c0, centered=FALSE)
    A.wt$strata <- as.numeric(A.wt$strata)
    for(trn in sort(unique(A.wt$strata))){
      fun <- stepfun(A.wt[A.wt$strata==trn,"time"], 
                     c(0,A.wt[A.wt$strata==trn,"hazard"]))
      A0$Haz[A0$Haz$trans==trn,"Haz"] <- fun(A0$Haz[A0$Haz$trans==trn,"time"])
    }
  }
  pt0 <- probtrans(A0, predt = s, variance=FALSE)[[h]][,c("time",
                                                            paste("pstate", j, sep=""))]
  return(pt0)
}


patp_b <- function(data, tmat, id2, h, j, s, times, wiaj_hat, n, weighted, LMAJ){
  if(LMAJ==FALSE){
    c0 <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data=data,
                method = "breslow")
    
    A0 <- msfit(object = c0, trans = tmat, variance=FALSE)
    if(weighted==TRUE){
      c0 <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), 
                  weights=(1/clust.size), data=data,
                  method = "breslow")
      A.wt <- basehaz(c0, centered=FALSE)
      A.wt$strata <- as.numeric(A.wt$strata)
      for(trn in sort(unique(A.wt$strata))){
        fun <- stepfun(A.wt[A.wt$strata==trn,"time"], 
                       c(0,A.wt[A.wt$strata==trn,"hazard"]))
        A0$Haz[A0$Haz$trans==trn,"Haz"] <- fun(A0$Haz[A0$Haz$trans==trn,"time"])
      }
    }
    
    P0 <- probtrans(A0, predt = s, 
                    variance=FALSE)[[h]][,c("time",
                                            paste("pstate", j, sep=""))]
  } else {
    P0 <- LMAJ2(msdata=data, tmat=tmat, id=id2, s=s, h=h, j=j, weighted=weighted)
  }

  P0_t <-stepfun(P0$time, c(P0[1,paste("pstate", j, sep="")], 
                            P0[,paste("pstate", j, sep="")]))                
  
  return(sqrt(n)*(P0_t(times)-wiaj_hat))
}



patp <- function(data, tmat, cid, id, h, j, s=0, weighted=FALSE,
                 LMAJ=FALSE, B=100, cband=FALSE){
  check.ic <- aggregate(data[,cid], by=list(data[,id]), 
                        FUN=sd, na.rm=TRUE)$x
  check.ic <- check.ic[!is.na(check.ic)]
  if(length(check.ic)>0){
    if(max(check.ic)>0){
      stop("Same unit(s) in more than 1 cluster (violation of the independent clusters assumption)")
    }
  }
  
  if(B<=0 & cband==TRUE){
    stop("Condidence bands cannot be caclulated based on <=0 bootstrap samples")
  } else if (B<1000 & cband==TRUE){
    warning("It is recommended to use at least 1000 bootstrap samples for confidence band calculation")
  }
  if(LMAJ==FALSE){
    c0 <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data=data,
                method = "breslow")
    
    A0 <- msfit(object = c0, trans = tmat, variance=FALSE)
    if(weighted==TRUE){
      ## msfit does not currently support weights and thus
      ## the weighted by cluster size cumulative transition
      ## intensities need to be manually inserted into A0
      M0 <- aggregate(rep(1,times=nrow(data)),
                      by = list( data[,cid], data[,id]),
                      FUN = mean)
      
      M <- aggregate(M0$x,
                     by = list( M0$Group.1),
                     FUN = sum)
      colnames(M) <- c(cid, "clust.size")
      data <- merge(data, M, by=cid)
      data <- data[order(data[,cid],data[,id]),]
      class(data) <- c("msdata", "data.frame")
      
      c0 <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), 
                  weights=(1/clust.size), data=data,
                  method = "breslow")
      A.wt <- basehaz(c0, centered=FALSE)
      A.wt$strata <- as.numeric(A.wt$strata)
      for(trn in sort(unique(A.wt$strata))){
        fun <- stepfun(A.wt[A.wt$strata==trn,"time"], 
                       c(0,A.wt[A.wt$strata==trn,"hazard"]))
        A0$Haz[A0$Haz$trans==trn,"Haz"] <- fun(A0$Haz[A0$Haz$trans==trn,"time"])
      }
    }
    P0 <- probtrans(A0, predt = s, 
                    variance=FALSE)[[h]][,c("time",
                                            paste("pstate", j, sep=""))]
  } else {
    if(weighted==TRUE){
      M0 <- aggregate(rep(1,times=nrow(data)),
                      by = list( data[,cid], data[,id]),
                      FUN = mean)
      
      M <- aggregate(M0$x,
                     by = list( M0$Group.1),
                     FUN = sum)
      colnames(M) <- c(cid, "clust.size")
      data <- merge(data, M, by=cid)
      data <- data[order(data[,cid],data[,id]),]
      class(data) <- c("msdata", "data.frame")
    }
    P0 <- LMAJ2(msdata=data, tmat=tmat, id=id, s=s, h=h, j=j, weighted=weighted)
  }

  colnames(P0) <- c("time", paste("P", h, j, sep=""))
  if(B==0){
    return(P0)
  } else {
    n <- length(unique(data[,cid]))
    boot <- msboot(patp_b, data=data, 
                   id=cid, B=B, verbose=0,
                   tmat=tmat, id2=id, h=h, j=j, s=s, times=P0$time,
                   wiaj_hat=P0[,paste("P", h, j, sep="")],
                   n=n, weighted=weighted, LMAJ=LMAJ)
    sigma <- apply(boot, 1, sd)
    se <- sigma/sqrt(n)
    
    ## cloglog transformation
    ll <- exp(-exp(log(-log(P0[,paste("P", h, j, sep="")]))-qnorm(0.975)*se/
                     (P0[,paste("P", h, j, sep="")]*
                        log(P0[,paste("P", h, j, sep="")]))))
    ul <- exp(-exp(log(-log(P0[,paste("P", h, j, sep="")]))+qnorm(0.975)*se/
                     (P0[,paste("P", h, j, sep="")]*
                        log(P0[,paste("P", h, j, sep="")]))))
    
    res <- cbind(P0, se, ll, ul)
    
    if(cband==TRUE){
      q_t <- 1/(1+sigma^2)
      
      jump.times <- P0$time[diff(c(P0[1,paste("P", h, j, sep="")],
                                   P0[,paste("P", h, j, sep="")]), lag=1)!=0]
      quant <- quantile(jump.times, probs=c(.05,.95))
      range <- (P0$time>=quant[1] & P0$time<=quant[2] & 
                  P0[,paste("P", h, j, sep="")]>0)
      
      B_t <- q_t[range]*boot[range,]/(log(P0[range, paste("P", h, j, sep="")])*
                                        P0[range, paste("P", h, j, sep="")])
      B_t <- abs(B_t)
      c_a <- apply(B_t, 2, max)
      c_a <- quantile(c_a, probs=0.95)
      
      ## cloglog transformation
      ll.band <- exp(-exp(log(-log(P0[,paste("P", h, j, sep="")]))+c_a/(sqrt(n)*q_t)))
      ll.band[!range] <- NA
      ul.band <- exp(-exp(log(-log(P0[,paste("P", h, j, sep="")]))-c_a/(sqrt(n)*q_t)))
      ul.band[!range] <- NA
      
      res <- cbind(res, ll.band, ul.band)
    }
    return(res)
  }
}



patp_test_b <- function(data, tmat, id2, h, j, s, group, times, D_hat, 
                        Wt, n, weighted, LMAJ){
  
  groups <- unique(data[,group])
  groups <- sort(groups[!is.na(groups)])
  
  Pt <- list()
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
      P0 <- probtrans(A0, predt = s, variance=FALSE)[[h]][,c("time",
                                                             paste("pstate", j, sep=""))]
    } else {
      P0 <- LMAJ2(msdata=data[data[,group]==g,], tmat=tmat, 
                  id=id2, s=s, h=h, j=j, weighted=weighted)
    }
    if(length(Pt)==0){
      Pt[[1]] <- stepfun(P0$time, c(P0[1,paste("pstate", j, sep="")], 
                                    P0[,paste("pstate", j, sep="")]))
    } else {
      Pt[[2]] <- stepfun(P0$time, c(P0[1,paste("pstate", j, sep="")], 
                                    P0[,paste("pstate", j, sep="")]))
    }
  }          
  D_boot <- Pt[[1]](times) - Pt[[2]](times)
  return(sqrt(n)*Wt*(D_boot-D_hat))
}



patp_test <- function(data, tmat, cid, id, group, h, j, s=0,
                      weighted=FALSE, LMAJ=FALSE, B=1000){
  check.ic <- aggregate(data[,cid], by=list(data[,id]), 
                        FUN=sd, na.rm=TRUE)$x
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
      
      P0 <- probtrans(A0, predt = s, variance=FALSE)[[h]][,c("time",
                                                             paste("pstate", j, sep=""))]
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


