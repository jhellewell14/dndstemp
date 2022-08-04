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

mod_test <- cmdstan_model("models/test.stan", 
                         cpp_options = list(stan_threads = TRUE))

fit <- mod_test$sample(
  data = list(pi_eq = rep(1/61, 61), mu = 1, omega = 1.27527, kappa = 2.72392, theta = 0.0878761,
              X = data_list$X, l = data_list$l),
  iter_warmup = 0, 
  iter_sampling = 1,
  threads_per_chain = 8,
  chains = 1, 
  fixed_param = TRUE
)

# Sum of likelihood over all positions
as.vector(fit$draws("lik"))
# Likelihood at each position
as.vector(fit$draws("likpos"))



