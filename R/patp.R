
## Function that calculates the working independence Aalen-Johansen estimator 
## of the population-averaged transition probabilities. Standard errors and
## and 95% confidence intervals and bands are also calculated.

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
