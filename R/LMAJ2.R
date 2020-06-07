
## Auxiliary function to create to calculate the landmark version
## of the working independence Aalen-Johansen estimator

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
