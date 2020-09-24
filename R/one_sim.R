## read in arguments from bash
in_args <- commandArgs(trailingOnly = TRUE)
data_id <- as.integer(in_args[1])

library(survival)
library(mstate)
source("LMAJ2.R")
source("patp_b.R")
source("patp_test_b.R")
source("patp_test.R")
source("patp.R")

mstat_long <- function(n = 100, a12 = 1, a13 = 0.5, a23 = 0.75, censor_rate=0.5)
{
  t12 <- rexp(n=n, rate=a12)
  t13 <- rexp(n=n, rate=a13)
  t23 <- rexp(n=n, rate=a23)
  c.t <- rexp(n=n, rate=censor_rate)
  cid <- as.integer(sort(sample(1:4, n, replace = TRUE)))
  # compare
  is_t12 <- t12 < c.t
  is_t13 <- t13 < c.t
  ## trans = 1
  is_ill <- t12 < t13 & t12 < c.t
  status <- ifelse(is_ill, 1, 0)
  tsub <- ifelse(is_t13, t13, c.t)
  Tstop <- ifelse(is_ill, t12, tsub)
  Tstart <- rep(0, n)
  ill_dat <- data.frame(
    id = seq(n), trans = 1, from = 1, to = 2,
    Tstart, Tstop, status, cid)
  ## trans = 2
  is_death <- t13 < t12 & t13 <c.t
  status <- ifelse(is_death, 1, 0)
  tsub <- ifelse(is_t12, t12, c.t)
  Tstop <- ifelse(is_death, t13, tsub)
  Tstart <- rep(0, n)
  death_dat <- data.frame(
    id = seq(n), trans = 2, from = 1, to = 3,
    Tstart, Tstop, status, cid)
  ## trans = 3
  is_23 <- t23 < c.t & is_ill
  status <- ifelse(is_23, 1, 0)
  Tstart <- t12
  Tstop <- ifelse(is_23, t12 + t23, c.t)
  Tstop[!status] <- NA
  d23_dat <-  data.frame(
    id = seq(n), trans = 3, from = 2, to = 3,
    Tstart, Tstop, status, cid)
  out <-  rbind(ill_dat, death_dat, d23_dat)
  out <- out[order(out$id),]
  out <- na.omit(out)
  out$time <- out$Tstop - out$Tstart
  class(out) <- c("msdata", class(out))
  return(out)
}


## True P12.
true_P12 <- function(a12, a13, a23, t){
  a12*(exp(-a23*t) - exp(-(a12+a13)*t))/(a12-a23+a13)
}

## Simulation function
sim_dat <- function(n = 100, tmat, a12 = 0.02, a13 = 0.05, a23 = 0.055,
                    censor_rate = 0.03) {
    mstat_long(n = n, a12 = a12, a13 = a13, a23 = a23,
               censor_rate = censor_rate)
}

run.sim <- function(dat_long, tmat,
                    a12 = 0.2, a13 = 0.5, a23 = 0.55,
                    times = c(0.5, 1, 2))
{
    P12_0 <- sapply(times, function(t) {
        true_P12(a12, a13, a23, t)
    })
    P12 <- patp(data = dat_long, tmat = tmat, cid = "cid", id = "id",
                h = 1, j = 2, s = 0, B = 1000)

    res <- sapply(times, function(t){P12[which.min(abs(P12$time - t)),]})
    P12_n = unlist(res[2, ])
    se_n = unlist(res[3,])
    ecp =  as.integer(P12_n - 1.96 * se_n <= P12_0 &
                      P12_n + 1.96 * se_n >= P12_0)
    res <- data.frame(P12_n, se_n, ecp )
    rownames(res) <- times
  return(res)
}

set.seed(123 + data_id)

## Define tmat
tmat <- transMat(x = list(c(2, 3), c(3), c()),
                 names = c("Health", "Illness" ,"Death"))
dat_long <- sim_dat(n = 200, tmat = tmat, a12 = 0.2, a13 = 0.5, a23 = 0.55,
                    censor_rate = 0.3)
out <- run.sim(dat_long, tmat = tmat, a12 = 0.2, a13 = 0.5, a23 = 0.55,
               times = c(0.5, 1, 1.5, 2))

## save the results
out_dir <- "output"
if (! dir.exists(out_dir)) {
  dir.create(out_dir)
}
saveRDS(out, file = file.path(out_dir, sprintf("one_sim_%d.rds", data_id)))
