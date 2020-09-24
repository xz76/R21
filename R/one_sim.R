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

sim_mstat <- function(n = 100, tmat, a12 = 0.02, a13 = 0.05, a23 = 0.055,
                      censor_rate = 0.03)
{
  t12 <- rexp(n=n, rate=a12)
  t13 <- rexp(n=n, rate=a13)
  t23 <- rexp(n=n, rate=a23)
  c.t <- rexp(n=n, rate=censor_rate)
  
  ill.s <- 1*(t12<=t13 & t12<=c.t)
  death <- ill.s*(t12+t23) + (1-ill.s)*t13
  
  death.s <- 1*(death <= c.t)
  ill <- ill.s*t12 + (1-ill.s)*ifelse(c.t<t13, c.t, t13)
  death <- death.s*death + (1 - death.s)*c.t
  
  id <- 1:n
  ## combine, sort and return
  out <- cbind(id, ill, ill.s, death, death.s)
  out <- as.data.frame(out)
  out$cid <- as.integer(sort(sample(1:4, n, replace = TRUE)))
  dat_long <- msprep(data = out, trans = tmat, time = c(NA, "ill", "death"),
                     status = c(NA, "ill.s", "death.s"), keep = c("cid", "id"))
  return(dat_long)
}

## True P12.
true_P12 <- function(a12, a13, a23, t){
  a12*(exp(-a23*t) - exp(-(a12+a13)*t))/(a12-a23+a13)
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
dat_long <- sim_mstat(n = 200, tmat = tmat, a12 = 0.2, a13 = 0.5, a23 = 0.55,
                    censor_rate = 0.3)
out <- run.sim(dat_long, tmat = tmat, a12 = 0.2, a13 = 0.5, a23 = 0.55,
               times = c(0.5, 1, 1.5, 2))

## save the results
out_dir <- "output"
if (! dir.exists(out_dir)) {
  dir.create(out_dir)
}
saveRDS(out, file = file.path(out_dir, sprintf("one_sim_%d.rds", data_id)))
