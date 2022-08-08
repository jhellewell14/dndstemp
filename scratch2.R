library(seqinr)
library(data.table)
library(cmdstanr)
source("generate_data.R")

data <- seqinr::read.alignment(file = "fits/porB/porB3.carriage.noindels.txt", format = "fasta")
data_list <- generate_data(data)
data_list$mu <- 1

# iteration loglikelihood loglik.seqs.     theta   kappa   omega0
# 1           5      -1957.71     -1958.48 0.1700000 1.00000 1.000000
# 3          15      -1939.37     -1940.16 0.1700000 1.63177 0.492084
# 151       755      -1934.01     -1934.22 0.1140420 2.83147 0.920608
dl <- list(pi_eq = rep(1/61, 61), mu = 1, omega = 1, kappa = 1, theta = 0.17,
           X = data_list$X, l = data_list$l, N = data_list$N, M = 42)
dl2 <- list(pi_eq = rep(1/61, 61), mu = 1, omega = 0.920608, kappa = 2.83147, theta = 0.1140420,
            X = data_list$X, l = data_list$l, N = data_list$N, M = 42)

nchains <- 1
cores <- (parallel::detectCores() - 2)
thr_per_chain <- floor(cores / nchains)
Sys.setenv(STAN_NUM_THREADS = 8)

mod_test <- cmdstan_model("models/optim.stan", 
                         cpp_options = list(stan_threads = TRUE))

fit <- mod_test$sample(
  data = dl2,
  iter_warmup = 500, 
  iter_sampling = 500,
  threads_per_chain = 8,
  chains = 1, 
  fixed_param = FALSE,
  refresh = 50
)

fit$summary(variables = c("omega", "kappa", "theta"))
as.vector(fit$draws("lik"))
sum(as.vector(fit$draws("v", format = "matrix")))
matrixStats::logSumExp(as.vector(fit$draws("v", format = "matrix")))
rbind(as.vector(fit$draws("v", format = "matrix"))[1:61], data_list$X[2,])
123 + 61
rbind(as.vector(fit$draws("v", format = "matrix"))[62:122], data_list$X[43,])

rbind(as.vector(fit$draws("v", format = "matrix"))[123:183], data_list$X[85,])
# Sum of likelihood over all positions
as.vector(fit$draws("lik"))[1]
