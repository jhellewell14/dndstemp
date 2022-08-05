library(seqinr)
library(data.table)
library(cmdstanr)
source("generate_data.R")

data <- seqinr::read.alignment(file = "fits/porB/porB3.carriage.noindels.txt", format = "fasta")
data_list <- generate_data(data)
data_list$mu <- 1

nchains <- 1
cores <- (parallel::detectCores() - 2)
thr_per_chain <- floor(cores / nchains)
Sys.setenv(STAN_NUM_THREADS = 8)
mod_mle <- cmdstan_model("models/mle.stan", 
                     cpp_options = list(stan_threads = TRUE))

# Find MLE
mod_mle$optimize(data = data_list, 
             algorithm = "lbfgs",
             init = initfn,
             threads = 8)

mod_test <- cmdstan_model("models/optim.stan", 
                         cpp_options = list(stan_threads = TRUE))

# 13.8

fit <- mod_test$sample(
  data = list(pi_eq = rep(1/61, 61), mu = 1, omega = 1, kappa = 1, theta = 0.17,
              X = data_list$X, l = data_list$l, N = data_list$N),
  iter_warmup = 0, 
  iter_sampling = 500,
  threads_per_chain = 8,
  chains = 1, 
  fixed_param = TRUE
)

# Sum of likelihood over all positions
as.vector(fit$draws("lik"))[1]
# Likelihood at each position
sum(as.vector(fit$draws("likposanc")))



