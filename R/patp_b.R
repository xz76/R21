
## Auxiliary function of nonparametric cluster bootstrapping.
## It is used for the calculation of standard errors and 
## 95% condidence intervals and bands.

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
