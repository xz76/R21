
## Auxiliary function for nonparametric cluster bootstrapping of the two-sample
## Kolmogorov-Smirnov-type test statistic.

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

