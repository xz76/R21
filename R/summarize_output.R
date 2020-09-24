get_point_est <- function(dat) {
    dat$P12_n
}
get_se_est <- function(dat) {
  dat$se_n
}
get_ecp_est <- function(dat){
  dat$ecp
}
dat <- readRDS("collected/one_sim.rds")
times = c(0.5, 1, 1.5, 2)
          
true_P12 <- function(a12, a13, a23, t){
  a12*(exp(-a23*t) - exp(-(a12+a13)*t))/(a12-a23+a13)
}
P12_0 <- sapply(times, function(t) {
  true_P12(a12 = 0.2, a13 = 0.5, a23 = 0.55, t)
})

est_mat <- sapply(dat, get_point_est)
bias <- rowMeans(est_mat) - P12_0
MCSD <- apply(est_mat, 1, sd)
ESE <- rowMeans(sapply(dat, get_se_est))
CP <- rowMeans(sapply(dat, get_ecp_est))

res <- as.data.frame(cbind(bias, MCSD, ESE, CP))
colnames(res) <- c("bias", "MCSD", "ESE", "CP")
rownames(res) <- times
