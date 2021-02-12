## Function to calculate the p-value for the two-sample comparison
## of transition probabilities based on a Kolmogorov-Smirnov-type
## test
library(matrixStats)
library(truncdist)
library(survival)
source("LMAJ2.R")
source("patp_test_b.R")
library("mstate")

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

      TstarTstop <- rweibull(1, shape=p12, scale=(1/(a12*exp(b.Z*Z))))
      Tstart3 <- rweibull(1, shape=p13, scale=(1/(a13*exp(b.Z*Z))))
      stop <- 1*((Tstart3 < TstarTstop) | (CC < TstarTstop)) #Stop if death or censoring first
      x1 <- min(c(TstarTstop, Tstart3)) # Exit time from state 1 (ignoring censoring)
      if (stop==0){
        Tstart <- 0
        Tstop <- x1
        from <- 1
        to <- 2
        while (stop==0){
          if(to[length(to)]==2){
            Tstop3 <- rtrunc(1, spec="weibull", a=Tstop[length(Tstop)], b=Inf,
                          shape=p23, scale=(1/(a23*exp(b.Z*Z))))
            Tstop1 <- rtrunc(1, spec="weibull", a=Tstop[length(Tstop)], b=Inf,
                          shape=p21, scale=(1/(a21*exp(b.Z*Z))))

            stop <- 1*((Tstop3 < Tstop1) | (CC < Tstop1)) #Stop if death or censoring first
            x2 <- min(c(Tstop1, Tstop3)) # Exit time from state 2 (ignoring censoring)

            Tstart <- c(Tstart, x1)
            Tstop <- c(Tstop, min(x2, CC))
            from <- c(from, 2)
            to <- c(to, ifelse(x2 <= CC, ifelse(Tstop1 <= Tstop3, 1, 3), 2))
          } else {
            TstarTstop <- rtrunc(1, spec="weibull", a=Tstop[length(Tstop)], b=Inf,
                          shape=p12, scale=(1/(a12*exp(b.Z*Z))))
            Tstart3 <- rtrunc(1, spec="weibull", a=Tstop[length(Tstop)], b=Inf,
                          shape=p13, scale=(1/(a13*exp(b.Z*Z))))

            stop <- 1*((Tstart3 < TstarTstop) | (CC < TstarTstop)) #Stop if death or censoring first
            x1 <- min(c(TstarTstop, Tstart3)) # Exit time from state 1 (ignoring censoring)

            Tstart <- c(Tstart, x2)
            Tstop <- c(Tstop, min(x1, CC))
            from <- c(from, 1)
            to <- c(to, ifelse(x1 <= CC, ifelse(TstarTstop <= Tstart3, 2, 3), 1))
          }
        }
      } else {
        Tstart <- 0
        Tstop <- min(x1, CC)
        from <- 1
        to <- ifelse(x1 < CC, 3, 1)
      }
    } else {
      Tstart <- 0
      Tstop <- NA
      from <- 0
      to <- NA
    }

    id <- rep(i, length(Tstart))
    Z <- rep(Z, length(Tstart))
    R <- rep(R, length(Tstart))
    outdata <- rbind(outdata, data.frame(id, Tstart, Tstop, from, to, Z, R))
  }
  # Tstart: start time of the interval
  # Tstop: end time of the interval
  # from: state at the start time of the interval
  # to: state at the end time of the interval
  # if from==to then right censoring occured at Tstop at this interval
  outdata <- outdata[order(outdata$id, outdata$Tstart, outdata$Tstop),]
  return(outdata)
}

#
# tmat <- transMat(x = list(c(2, 3), c(1, 3), c()),
#                  names = c("Health", "Illness", "Death"))
reshape_long <- function(dat, tmat) {
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
    is_censor <- sub_dat$from == sub_dat$to

    out <- lapply(seq_along(is_censor), function(i) {
      this_censor <- is_censor[i]
      if (this_censor) {
        res <- data.frame(id = id,
                          Tstart = sub_dat$Tstart[i],
                          Tstop = sub_dat$Tstop[i],
                          time = sub_dat$Tstop[i] - sub_dat$Tstart[i],
                          from = sub_dat$from[i])
        res <- merge(res, tmp_dat, by = "from")
        res <- res[, c("id", "from", "to", "Tstart", "Tstop", "time",
                       "trans")]
        res$status <- 0
      } else {
        res <- data.frame(id = id,
                          Tstart = sub_dat$Tstart[i],
                          Tstop = sub_dat$Tstop[i],
                          from = sub_dat$from[i],
                          to = sub_dat$to[i],
                          time = sub_dat$Tstop[i] - sub_dat$Tstart[i],
                          status = 1)
        res <- merge(res, tmp_dat, by = c("from", "to"))
        res <- res[, c("id", "from", "to", "trans", "Tstart",
                       "Tstop", "time", "status")]
      }
      res
    })
    do.call(rbind, out)
  }
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
sop_t <- function(data, tau=NULL, S, T_c, ipw=0, trans,
                     times= NULL){
  if(is.null(tau)){
    tau <- max(data$Tstop)
  }
  if(is.null(times)){

  }
  CTI <- list()
  counter <- 1
  for(h in T_c){
    for(j in S[trans[h,]]){
      data_h <- data[data$from==h,]
      data_h$delta <- 1*(data_h$to==j)
      if(ipw==0){
        fit <- coxph(Surv(Tstart,Tstop,delta,type="counting")~1, data=data_h)
      } else {
        fit <- coxph(Surv(Tstart,Tstop,delta,type="counting")~1,
                     weight=weightVL, data=data_h)
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

  tt<-sort(unique(data[data$to != data$from,"Tstop"]))
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
  p_0 <- sapply(S, function(i){sum(data[data$Tstart==0,"from"] == i)/nrow(data[data$Tstart==0,])})
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

D.AUC <- function(data,cov,tau,S,T_c,ipw, trans){
  data <- tmp
  cov <- data$group
  tau = NULL
  T_c = 1:2
  S = 1:3
  ipw = 0
  if (cov != 0 && cov != 1) {
    stop("The 'cov' has to be 0 or 1.")
  }
  # P1 <- sop_t(data[data[,cov]==1,], tau=tau, S=S, T_c=T_c, ipw=ipw, trans = tmatrix)
  # P0 <- sop_t(data[data[,cov]==0,], tau=tau, S=S, T_c=T_c, ipw=ipw, trans = tmatrix)

  P1 <- sop_t(data[cov==1,], tau=tau, S=S, T_c=T_c, ipw=ipw, trans = tmatrix)
  P0 <- sop_t(data[cov==0,], tau=tau, S=S, T_c=T_c, ipw=ipw, trans = tmatrix)

  dAUC <- vector(length = 5)
  for(j in S){
    p1 <- P1[,c("t",paste("p",j,sep=""))]
    p0 <- P0[,c("t",paste("p",j,sep=""))]

    jump1 <- (p1[,paste("p",j,sep="")]!=c(0,p1[1:(nrow(p1)-1),paste("p",j,sep="")]))
    jump0 <- (p0[,paste("p",j,sep="")]!=c(0,p0[1:(nrow(p0)-1),paste("p",j,sep="")]))

    tt<-sort(unique(c(p1[jump1,"t"],p0[jump0,"t"])))

    dm<-diff(c(tt,tau),lag=1)


    elemenTstart<-function(t){
      max(1:length(p1$t)*(p1$t<=t))
    }
    elementfrom<-sapply(tt,elemenTstart)
    elementfrom<-(elementfrom==0)+(elementfrom>0)*elementfrom

    element0<-function(t){
      max(1:length(p0$t)*(p0$t<=t))
    }
    elements0<-sapply(tt,element0)
    elements0<-(elements0==0)+(elements0>0)*elements0

    D_t <- p1[elementfrom,paste("p",j,sep="")]-p0[elements0,paste("p",j,sep="")]
    dAUC[j] <- sum(D_t*dm)
  }
  dAUC <- c(dAUC, (dAUC[1] + dAUC[2] + dAUC[3]))

  p1 <- P1[,"p3"]/(rowSums(P1[,c("p1","p2","p3")]))
  p1 <- cbind(P1$t,p1)
  p1 <- as.data.frame(p1)
  colnames(p1) <- c("t","p3")
  p0 <- P0[,"p3"]/(rowSums(P0[,c("p1","p2","p3")]))
  p0 <- cbind(P0$t,p0)
  p0 <- as.data.frame(p0)
  colnames(p0) <- c("t","p3")

  jump1 <- (p1[,"p3"]!=c(0,p1[1:(nrow(p1)-1),"p3"]))
  jump0 <- (p0[,"p3"]!=c(0,p0[1:(nrow(p0)-1),"p3"]))

  tt<-sort(unique(c(p1[jump1,"t"],p0[jump0,"t"])))

  dm<-diff(c(tt,tau),lag=1)

  elemenTstart<-function(t){
    max(1:length(p1$t)*(p1$t<=t))
  }
  elementfrom<-sapply(tt,elemenTstart)
  elementfrom<-(elementfrom==0)+(elementfrom>0)*elementfrom

  element0<-function(t){
    max(1:length(p0$t)*(p0$t<=t))
  }
  elements0<-sapply(tt,element0)
  elements0<-(elements0==0)+(elements0>0)*elements0

  D_t <- p1[elementfrom,"p3"]-p0[elements0,"p3"]
  dAUC <- c(dAUC, sum(D_t*dm))
  return(dAUC)
}

patp_test <- function(data, tmat,cid = "cid", id = "id", group = "group", h = 1,
                      j = 2, s=0, weighted=FALSE, LMAJ=FALSE, B=1000, ipw = 0){
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
    if (!LMAJ) {
      # c0 <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans),
      #             data=data[data[,group]==g,], method = "breslow")
      #
      # A0 <- msfit(object = c0, trans = tmat, variance=FALSE)
      # if(weighted==TRUE){
      #   c0 <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans),
      #               weights=(1/clust.size), data=data[data[,group]==g,],
      #               method = "breslow")
      #   A.wt <- basehaz(c0, centered=FALSE)
      #   A.wt$strata <- as.numeric(A.wt$strata)
      #   for(trn in sort(unique(A.wt$strata))){
      #     fun <- stepfun(A.wt[A.wt$strata==trn,"time"],
      #                    c(0,A.wt[A.wt$strata==trn,"hazard"]))
      #     A0$Haz[A0$Haz$trans==trn,"Haz"] <- fun(A0$Haz[A0$Haz$trans==trn,"time"])
      #   }
      # }
     # P0 <- probtrans(A0, predt = s,
     #                 variance=FALSE)[[h]][,c("time", paste("pstate", j, sep=""))]
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

    } else {
      Pt[[2]] <- stepfun(P0$t, rep1(P0[, paste0("p", j)]))
      tms[[2]] <- P0$t[diff(rep1(P0[, paste0("p", j)]), lag=1)!=0]
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

  Wt <- rowProds(EY)/rowSums(EY)
 # Wt <- rowprods(EY)/rowSums(EY)
  tms <- tms[!is.na(Wt)]
  if(length(tms)==0 | max(Wt, na.rm=TRUE)==0){
    stop("Weights NA or 0 for all timepoints")
  }
  Wt <- Wt[!is.na(Wt)]
  D_hat=(Pt[[1]](tms) - Pt[[2]](tms))

  ### KS result p-value = 1 ####

  patp_test_b <- function(data, tmat, h, j, group, times, D_hat,
                          Wt, n, weighted, LMAJ = FALSE){
    groups <- unique(data[,group])
    groups <- sort(groups[!is.na(groups)])
    Pt <- lapply(seq_along(groups), function(k) {
      g <- groups[k]
      if (! LMAJ) {
        group_data <- data[data[,group]== g,]
        P0 <- sop_t(data = group_data, tau=NULL, S = seq(nrow(tmat)), T_c = seq(nrow(tmat) - 1),
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

  diff_boot <- function ( data, B, id,  verbose = 0, 
                          tmat, group, h=h, j=j,
                          times, D_hat, Wt,
                          n, weighted, LMAJ=FALSE) 
  { 
    ids <- unique(data[[id]])
    n <- length(ids)
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
  
  Diff_boot <- diff_boot(data=data, id = cid,  B=B, verbose=0,
                         tmat=tmat, group=group, h=h, j=j,
                         times=tms, D_hat=D_hat, Wt=Wt,
                         n=n, weighted=weighted, LMAJ=LMAJ)

  KS <- max(abs(sqrt(n)*Wt*D_hat))
  KS.b <- apply(abs(Diff_boot),2,max)
  pval <- mean(KS.b >= KS)
  names(pval) <- "p-value"
  return(pval)
}

## Example
tmatrix <- trans(state_names = c("health", "illness", "death"),from = c(1, 1, 1, 2, 2),
                 to = c(2, 2, 3, 3, 1))
set.seed(0212)
tmp <- simulate(300)
tmp <- reshape_long(tmp, tmatrix)
tmp <- sim_group(tmp, cid = 5)
KS_res <- patp_test(tmp, tmatrix, cid = "cid", id = "id", group = "group", h = 1, j = 2, B = 1000)
