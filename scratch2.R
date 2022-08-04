library(seqinr)
library(data.table)
library(cmdstanr)
source("generate_data.R")


data <- seqinr::read.alignment(file = "fits/porB/porB3.carriage.noindels.txt", format = "fasta")
data_list <- generate_data(data)

mod <- cmdstan_model("models/test.stan", 
                     cpp_options = list(stan_threads = TRUE))


fit <- mod$sample(
  data = list(pi_eq = rep(1/61, 61), mu = 1, omega = 1.27527, kappa = 2.72392, theta = 0.0878761,
              X = data_list$X, l = data_list$l),
  iter_warmup = 0, 
  iter_sampling = 1,
  threads_per_chain = 8,
  chains = 1, 
  fixed_param = TRUE
)
y <- as.vector(fit$draws("likpos"))

as.vector(fit$draws("lik"))
# Likelihood at each position

as.vector(fit$draws("phi"))
as.vector(fit$draws("likposanc"))
# Sum of likelihood over all positions


m_AB <- matrix(fit$draws(variables = "m_AB", format = "matrix"), byrow = TRUE, nrow = 61)
mutmat <- matrix(fit$draws(variables = "mutmat", format = "matrix"), byrow = TRUE, nrow = 61)
V <- matrix(fit$draws(variables = "V", format = "matrix"), byrow = TRUE, nrow = 61)
Ve <- as.vector(fit$draws(variables = "Ve", format = "matrix"))
# plot(x = 1:61,y = rowMeans(replicate(expr = gtools::rdirichlet(n = 1, alpha = m_AB[1,]),                                      
#                                      n = 1000, simplify = "matrix")), ylim = c(0,1))
# fit$draws("log_lik")
# fit$summary()
# 
# mr <- as.vector(fit$draws(variables = "meanrate", format = "matrix"))
# D <- as.vector(fit$draws(variables = "D", format = "matrix"))
# E <- as.vector(fit$draws(variables = "E", format = "matrix"))
# M <- matrix(as.vector(fit$draws(variables = "mutmat", format = "matrix")), byrow = TRUE, nrow = 61)

# logvec <- as.vector(fit$draws(variables = "logvec", format = "matrix"))
# 
# plot(x = 1:61,y = rowMeans(replicate(expr = gtools::rdirichlet(n = 1, alpha = m_AB[57,]),                                      
#                                      n = 1000, simplify = "matrix")), ylim = c(0,1))
# points(1:61, data_list$X[lc,]/sum(data_list$X[lc,]), col = "red")
# 
# rbind(m_AB[57,], data_list$X[lc,])
# 
fitfn <- function(x){
  print(x)
  fit <- mod$sample(
    data = list(theta = x[1], omega = x[2], kappa = x[3], pi_eq = rep(1/61, 61), mu = 1,
                X = data_list$X, l = data_list$l),
    iter_warmup = 0,
    iter_sampling = 1,
    threads_per_chain = 1,
    chains = 1,
    fixed_param = TRUE
  )
  return(as.vector(fit$draws("lik", format = "matrix")))
}


optim(par = c(1, 1, 1), fn = fitfn,
      control = list(fnscale = -1),
      lower = c(0.01, 0.1, 0.1),
      upper = c(2, 5, 5),
      method = "L-BFGS-B")
# 
# ?optim
# M <- matrix(as.vector(fit$draws(variables = "mutmat", format = "matrix")), byrow = TRUE, nrow = 61)
# # fit$draws("meanrate")
# # V <- matrix(as.vector(fit$draws(variables = "V", format = "matrix")), byrow = TRUE, nrow = 61)
# # Va <- matrix(as.vector(fit$draws(variables = "Va", format = "matrix")), byrow = TRUE, nrow = 61)
# # V_inv <- matrix(as.vector(fit$draws(variables = "V_inv", format = "matrix")), byrow = TRUE, nrow = 61)
# # V_inv <- matrix(as.vector(fit$draws(variables = "V_inv", format = "matrix")), byrow = TRUE, nrow = 61)
# # D <- as.vector(fit$draws(variables = "D", format = "matrix"))
# # Rt <- matrix(fit$draws(variables = "Rt", format = "matrix"), byrow = TRUE, nrow = 61)
# 
# 
# dmlpmf(y = data_list$X[lc,], alpha = as.vector(fit$draws(variables = "m_AB", format = "matrix")))
# 
# 
# fit$draws(variables = "mutmat", format = "matrix")
# M1 <- matrix(as.vector(fit$draws(variables = "mutmat", format = "matrix")), byrow = FALSE, nrow = 61)
# M2 <- PDRM(mu = 9, omega = 3, kappa = 2, pi_eq = rep(1, 61))
# all.equal(M1, M2)
# 
# stl <- function(mu, omega, kappa, an, lc){
#   fit <- mod$sample(
#     data = list(mu = mu, omega = omega, kappa = kappa, pi_eq = rep(1/61, 61),
#                 k = an, loc = lc, X = data_list$X, l = data_list$l),
#     iter_warmup = 0, 
#     iter_sampling = 1,
#     threads_per_chain = 1,
#     chains = 1, 
#     fixed_param = TRUE
#   )
#   return(as.vector(fit$draws(variables = "logvec", format = "matrix")))
# }
# c(7, 23, 39, 52)
# stl(mu = 3, omega = 1, kappa = 1, lc = 10, an = 8)
# stl(mu = 3, omega = 1, kappa = 1, lc = 10, an = 23)
# stl(mu = 3, omega = 1, kappa = 1, lc = 10, an = 39)
# stl(mu = 3, omega = 1, kappa = 1, lc = 10, an = 52)
# 
# M2[1:5, 1:5]
# M_col <- matrix("0", 61, 61)
# M_col[M2 == M2[1,2]] <- "S"
# M_col[M2 == M2[1,3]] <- "NS"
# M_col[M2 == M2[1,5]] <- "NS"
# M_col[M2 == M2[7,5]] <- "S"
# M_col[M2 < 0] <- "D"
# M_col[M2 == 0] <- "o"
# 
# plot(M_col, col = c("purple","red", "white", "blue"))
# 
# colnames(M_col) <- rev(tripletNames_noSTO)
# rownames(M_col) <- rev(tripletNames_noSTO)
# 
# M_rev <- matrix(rev(M_col), 61, 61)
# plot(M_rev, col = c("purple","red", "white", "blue"))
# 
# PDRM <- function(mu, kappa, omega, pi_eq) {
#   M = matrix(0, nrow = 61, ncol = 61);
#   
#   M[1,2] = kappa
#   M[1,3] = omega
#   M[1,4] = omega
#   M[1,5] = kappa*omega
#   M[1,9] = omega
#   M[1,11] = omega
#   M[1,14] = kappa*omega
#   M[1,30] = omega
#   M[1,46] = omega
#   M[2,3] = omega
#   M[2,4] = omega
#   M[2,6] = kappa*omega
#   M[2,10] = omega
#   M[2,12] = omega
#   M[2,15] = kappa*omega
#   M[2,31] = omega
#   M[2,47] = omega
#   M[3,4] = kappa
#   M[3,7] = kappa*omega
#   M[3,16] = kappa
#   M[3,32] = omega
#   M[3,48] = omega
#   M[4,8] = kappa*omega
#   M[4,13] = omega
#   M[4,17] = kappa
#   M[4,33] = omega
#   M[4,49] = omega
#   M[5,6] = kappa
#   M[5,7] = 1
#   M[5,8] = 1
#   M[5,9] = omega
#   M[5,11] = omega
#   M[5,18] = kappa*omega
#   M[5,34] = omega
#   M[5,50] = omega
#   M[6,7] = 1
#   M[6,8] = 1
#   M[6,10] = omega
#   M[6,12] = omega
#   M[6,19] = kappa*omega
#   M[6,35] = omega
#   M[6,51] = omega
#   M[7,8] = kappa
#   M[7,20] = kappa*omega
#   M[7,36] = omega
#   M[7,52] = omega
#   M[8,13] = omega
#   M[8,21] = kappa*omega
#   M[8,37] = omega
#   M[8,53] = omega
#   M[9,10] = kappa
#   M[9,11] = kappa*omega
#   M[9,22] = kappa*omega
#   M[9,38] = omega
#   M[9,54] = omega
#   M[10,12] = kappa*omega
#   M[10,23] = kappa*omega
#   M[10,39] = omega
#   M[10,55] = omega
#   M[11,12] = kappa
#   M[11,13] = omega
#   M[11,26] = kappa*omega
#   M[11,42] = omega
#   M[11,58] = omega
#   M[12,13] = omega
#   M[12,27] = kappa*omega
#   M[12,43] = omega
#   M[12,59] = omega
#   M[13,29] = kappa*omega
#   M[13,45] = omega
#   M[13,61] = omega
#   M[14,15] = kappa
#   M[14,16] = 1
#   M[14,17] = 1
#   M[14,18] = kappa*omega
#   M[14,22] = omega
#   M[14,26] = omega
#   M[14,30] = omega
#   M[14,46] = omega
#   M[15,16] = 1
#   M[15,17] = 1
#   M[15,19] = kappa*omega
#   M[15,23] = omega
#   M[15,27] = omega
#   M[15,31] = omega
#   M[15,47] = omega
#   M[16,17] = kappa
#   M[16,20] = kappa*omega
#   M[16,24] = omega
#   M[16,28] = omega
#   M[16,32] = omega
#   M[16,48] = omega
#   M[17,21] = kappa*omega
#   M[17,25] = omega
#   M[17,29] = omega
#   M[17,33] = omega
#   M[17,49] = omega
#   M[18,19] = kappa
#   M[18,20] = 1
#   M[18,21] = 1
#   M[18,22] = omega
#   M[18,26] = omega
#   M[18,34] = omega
#   M[18,50] = omega
#   M[19,20] = 1
#   M[19,21] = 1
#   M[19,23] = omega
#   M[19,27] = omega
#   M[19,35] = omega
#   M[19,51] = omega
#   M[20,21] = kappa
#   M[20,24] = omega
#   M[20,28] = omega
#   M[20,36] = omega
#   M[20,52] = omega
#   M[21,25] = omega
#   M[21,29] = omega
#   M[21,37] = omega
#   M[21,53] = omega
#   M[22,23] = kappa
#   M[22,24] = omega
#   M[22,25] = omega
#   M[22,26] = kappa*omega
#   M[22,38] = omega
#   M[22,54] = omega
#   M[23,24] = omega
#   M[23,25] = omega
#   M[23,27] = kappa*omega
#   M[23,39] = omega
#   M[23,55] = omega
#   M[24,25] = kappa
#   M[24,28] = kappa*omega
#   M[24,40] = omega
#   M[24,56] = omega
#   M[25,29] = kappa*omega
#   M[25,41] = omega
#   M[25,57] = omega
#   M[26,27] = kappa
#   M[26,28] = 1
#   M[26,29] = 1
#   M[26,42] = omega
#   M[26,58] = omega
#   M[27,28] = 1
#   M[27,29] = 1
#   M[27,43] = omega
#   M[27,59] = omega
#   M[28,29] = kappa
#   M[28,44] = 1
#   M[28,60] = omega
#   M[29,45] = 1
#   M[29,61] = omega
#   M[30,31] = kappa
#   M[30,32] = 1
#   M[30,33] = omega
#   M[30,34] = kappa*omega
#   M[30,38] = omega
#   M[30,42] = omega
#   M[30,46] = kappa*omega
#   M[31,32] = 1
#   M[31,33] = omega
#   M[31,35] = kappa*omega
#   M[31,39] = omega
#   M[31,43] = omega
#   M[31,47] = kappa*omega
#   M[32,33] = kappa*omega
#   M[32,36] = kappa*omega
#   M[32,40] = omega
#   M[32,44] = omega
#   M[32,48] = kappa*omega
#   M[33,37] = kappa*omega
#   M[33,41] = omega
#   M[33,45] = omega
#   M[33,49] = kappa*omega
#   M[34,35] = kappa
#   M[34,36] = 1
#   M[34,37] = 1
#   M[34,38] = omega
#   M[34,42] = omega
#   M[34,50] = kappa*omega
#   M[35,36] = 1
#   M[35,37] = 1
#   M[35,39] = omega
#   M[35,43] = omega
#   M[35,51] = kappa*omega
#   M[36,37] = kappa
#   M[36,40] = omega
#   M[36,44] = omega
#   M[36,52] = kappa*omega
#   M[37,41] = omega
#   M[37,45] = omega
#   M[37,53] = kappa*omega
#   M[38,39] = kappa
#   M[38,40] = omega
#   M[38,41] = omega
#   M[38,42] = kappa*omega
#   M[38,54] = kappa*omega
#   M[39,40] = omega
#   M[39,41] = omega
#   M[39,43] = kappa*omega
#   M[39,55] = kappa*omega
#   M[40,41] = kappa
#   M[40,44] = kappa*omega
#   M[40,56] = kappa*omega
#   M[41,45] = kappa*omega
#   M[41,57] = kappa*omega
#   M[42,43] = kappa
#   M[42,44] = omega
#   M[42,45] = omega
#   M[42,58] = kappa*omega
#   M[43,44] = omega
#   M[43,45] = omega
#   M[43,59] = kappa*omega
#   M[44,45] = kappa
#   M[44,60] = kappa*omega
#   M[45,61] = kappa*omega
#   M[46,47] = kappa
#   M[46,48] = 1
#   M[46,49] = 1
#   M[46,50] = kappa*omega
#   M[46,54] = omega
#   M[46,58] = omega
#   M[47,48] = 1
#   M[47,49] = 1
#   M[47,51] = kappa*omega
#   M[47,55] = omega
#   M[47,59] = omega
#   M[48,49] = kappa
#   M[48,52] = kappa*omega
#   M[48,56] = omega
#   M[48,60] = omega
#   M[49,53] = kappa*omega
#   M[49,57] = omega
#   M[49,61] = omega
#   M[50,51] = kappa
#   M[50,52] = 1
#   M[50,53] = 1
#   M[50,54] = omega
#   M[50,58] = omega
#   M[51,52] = 1
#   M[51,53] = 1
#   M[51,55] = omega
#   M[51,59] = omega
#   M[52,53] = kappa
#   M[52,56] = omega
#   M[52,60] = omega
#   M[53,57] = omega
#   M[53,61] = omega
#   M[54,55] = kappa
#   M[54,56] = omega
#   M[54,57] = omega
#   M[54,58] = kappa*omega
#   M[55,56] = omega
#   M[55,57] = omega
#   M[55,59] = kappa*omega
#   M[56,57] = kappa
#   M[56,60] = kappa*omega
#   M[57,61] = kappa*omega
#   M[58,59] = kappa
#   M[58,60] = 1
#   M[58,61] = 1
#   M[59,60] = 1
#   M[59,61] = 1
#   M[60,61] = kappa
#   
#   M <- M + t(M)
#   
#   M <- mu * M
#   
#   for (i in 1:61) {
#     M[i,] <- M[i,] * pi_eq
#   }
#   
#   for (i in 1:61){
#     M[i, i] = -sum(M[i,])
#   }
#   
#   return(M)
# }
# 
# dmlpmf <-function(y, alpha) {
#   alpha_plus <- sum(alpha)
#   return(lgamma(alpha_plus) + sum(lgamma(alpha + y))
#          - lgamma(alpha_plus+sum(y)) - sum(lgamma(alpha)))
# }
# 
# 
# 
# 
# 
# mutmat <- PDRM(mu = 1,kappa =  2.5,omega = 1, pi_eq = rep(1/61, 61))
# which(data_list$X[1, ] > 0)
# sum(llcalc(mu = 4,kappa =  3,omega = 1, l = 1))

