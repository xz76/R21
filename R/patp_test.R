patp_test <- function(data, tmat, cid = "cid", id = "id", group = "group", h = 1,
                      j = 2, s = 0, weighted=FALSE, LMAJ=FALSE, B = 1000, ipw = 0,
                      method = "linear", tau = NULL){
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
    with(data, tapply(rep(1,times=nrow(data)),
                            paste(cid, group, id, sep = ":"), FUN = mean))
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
  tmax <- list()
  S = seq(nrow(tmat))
  T_c = seq(nrow(tmat) - 1)
  for(g in groups){
    if (!LMAJ) {
      group_data <- data[data[,group]==g,]
      P0 <- sop_t(data = group_data, tau=NULL, S = seq(nrow(tmat)), T_c = seq(nrow(tmat) - 1),
                  ipw=0, trans = tmat, times = NULL)
      
    } else {
      P0 <- LMAJ2(msdata=data[data[,group]==g,], tmat=tmat,
                  id=id, s=s, h=h, j=j, weighted=weighted)
    }
    rep1 <- function(x){
      x[c(1, seq_along(x))]
    }
    if(length(Pt)==0){
      Pt[[1]] <- stepfun(P0$t, rep1(P0[, paste0("p", j)]))
      tms[[1]] <- P0$t[diff(rep1(P0[, paste0("p", j)]), lag=1)!=0]
      tmax[[1]] <- max(P0[P0[, paste0("p",j)]>0, "t"])
    } else {
      Pt[[2]] <- stepfun(P0$t, rep1(P0[, paste0("p", j)]))
      tms[[2]] <- P0$t[diff(rep1(P0[, paste0("p", j)]), lag=1)!=0]
      tmax[[2]] <- max(P0[P0[, paste0("p",j)]>0, "t"])
    }
  }
  
  tms <- sort(unique(c(tms[[1]], tms[[2]])))
  tau0 <- max(data$Tstop)
  tmax <- min(tmax[[1]], tmax[[2]])
  tms <- tms[tms <= tmax]
  tau <- min(tau0, tmax)
  
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
  
  Wt <- rowProds(EY)/rowSums(EY)
  
  tms <- tms[!is.na(Wt)]
  if(length(tms)==0 | max(Wt, na.rm=TRUE)==0){
    stop("Weights NA or 0 for all timepoints")
  }
  Wt <- Wt[!is.na(Wt)]
  D_hat=(Pt[[1]](tms) - Pt[[2]](tms))
  
  estimator <- function(data,cov,tau,S,T_c,ipw, trans, Wt, method){
    cov <- data$group
    ipw <- 0
    if (cov != 0 && cov != 1) {
      stop("The 'cov' has to be 0 or 1.")
    }
    
    P1 <- sop_t(data[cov==1,], tau=tau, S=S, T_c=T_c, ipw=ipw, trans = tmat)
    P0 <- sop_t(data[cov==0,], tau=tau, S=S, T_c=T_c, ipw=ipw, trans = tmat)
    res <- vector(length = length(S))
    for(j in S){
      p1 <- P1[,c("t",paste("p",j,sep=""))]
      p0 <- P0[,c("t",paste("p",j,sep=""))]
      
      dm<-diff(c(tms,tau),lag=1)
      elemenTstart<-function(t){
        max(1:length(p1$t)*(p1$t<=t))
      }
      elementfrom<-sapply(tms,elemenTstart)
      elementfrom<-(elementfrom==0)+(elementfrom>0)*elementfrom
      
      element0<-function(t){
        max(1:length(p0$t)*(p0$t<=t))
      }
      elements0<-sapply(tms,element0)
      elements0<-(elements0==0)+(elements0>0)*elements0
      D_t <- p1[elementfrom,paste("p",j,sep="")]-p0[elements0,paste("p",j,sep="")]
      if (method == "KS"){
        res[j] <- max(abs(Wt * D_t * sqrt(n)))
      }
      if (method == "linear"){
        res[j] <- sum(Wt * D_t * dm)
      }
    }
    return(res)
  }
  
  if (method == "KS") {
    diff_boot <- function (data, B, id, verbose = 0, tmat, group, h=h, j=j, 
                           times, D_hat, Wt, n, weighted, LMAJ)
    {
      ids <- unique(data[[id]])
      n <- length(ids)
      patp_test_b <- function(data, tmat, h, j, group, times, D_hat,
                              Wt, n, weighted, LMAJ = FALSE){
        groups <- unique(data[,group])
        groups <- sort(groups[!is.na(groups)])
        Pt <- lapply(seq_along(groups), function(k) {
          g <- groups[k]
          if (! LMAJ) {
            group_data <- data[data[,group]== g,]
            P0 <- sop_t(data = group_data, tau=NULL, S = seq(nrow(tmat)), 
                        T_c = seq(nrow(tmat) - 1),
                        ipw=0, trans = tmat, times = NULL)
            
          } else {
            P0 <- LMAJ2(msdata=data[data[,group]==g,], tmat=tmat,
                        id=id, s=s, h=h, j=j, weighted=weighted)
          }
          stepfun(P0$t, rep1(P0[, paste0("p", j)]))
        })
        D_boot <- Pt[[1]](times) - Pt[[2]](times)
        return(sqrt(n)*Wt*(D_boot-D_hat))
      }
      th <- patp_test_b(data, tmat=tmat, group=group, h=h, j=j,
                        times=tms, D_hat=D_hat, Wt=Wt,
                        n=n, weighted=weighted, LMAJ=FALSE)
      res <- matrix(NA, length(th), B)
      for (b in 1:B) {
        if (verbose > 0) {
          cat("\nBootstrap replication", b, "\n")
          flush.console()
        }
        bootdata <- NULL
        bids <- sample(ids, replace = TRUE)
        bidxs <- unlist(sapply(bids, function(x) which(x == data[[id]])))
        bootdata <- data[bidxs, ]
        if (verbose > 0) {
          print(date())
          print(events(bootdata))
          cat("applying theta ...")
        }
        thstar <- patp_test_b(bootdata, tmat=tmat, group=group, h=h, j=j,
                              times=tms, D_hat=D_hat, Wt=Wt,
                              n=n, weighted=weighted, LMAJ=FALSE)
        res[, b] <- thstar
      }
      if (verbose)
        cat("\n")
      return(res)
    }
    
    Diff_boot <- diff_boot(data = data, id = cid,  B = B, verbose = 0,
                           tmat = tmat, group = group, h = h, j = j,
                           times = tms, D_hat = D_hat, Wt = Wt,
                           n = n, weighted = weighted, LMAJ = LMAJ)
    KS <- estimator(data = data,cov = data$group ,tau = NULL ,S = S,
                    T_c = T_c, ipw = 0, trans = tmat, Wt = Wt, method = "KS")
    KS.b <- apply(abs(Diff_boot),2,max)
    pval <- mean(KS.b >= KS[j])
    names(pval) <- "p-value"
    return(pval)
  }
  
  if (method == "linear") {
    
    Z0 <-  estimator( data = data, cov = data$group, tau = NULL, S = S,
                      T_c = T_c, ipw = 0, trans = tmat, Wt = Wt, method = "linear")
    dauc_boot <- function (data, cov, tau =NULL, S, T_c, ipw, trans,
                           Wt, B,  verbose = 0, id)
    {
      ids <- unique(data[[id]])
      n <- length(ids)
      th <- estimator( data = data, cov = data$group, tau = tau, S = S,
                       T_c = T_c, ipw = 0, trans = tmat, Wt = Wt, method = "linear")
      res <- matrix(NA, length(th), B)
      for (b in 1:B) {
        if (verbose > 0) {
          cat("\nBootstrap replication", b, "\n")
          flush.console()
        }
        bootdata <- NULL
        bids <- sample(ids, replace = TRUE)
        bidxs <- unlist(sapply(bids, function(x) which(x == data[[id]])))
        bootdata <- data[bidxs, ]
        if (verbose > 0) {
          print(date())
          print(events(bootdata))
          cat("applying theta ...")
        }
        thstar <-  estimator( data = bootdata, cov = bootdata$group, tau = tau, S = S,
                              T_c = T_c, ipw = ipw, trans = tmat, Wt = Wt, method = "linear")
        res[, b] <- thstar
      }
      if (verbose)
        cat("\n")
      return(res)
    }
    dauc_res <- dauc_boot( data = data, cov = data$group, S = S, T_c = T_c, 
                           ipw = ipw, trans = tmat, Wt = Wt, B = B, id = cid)
    dauc_sd <- apply(dauc_res, 1, sd)
    
    T_linear <- Z0/dauc_sd
    
    pval_lin <- 2*pnorm(abs(T_linear), lower.tail = FALSE)
    names(pval_lin) <- "p-value"
    return(pval_lin)
  }
}