library(cmdstanr)
source("generate_data.R")
data_list <- readRDS("data_list.RDS")

cores <- 4
nchains <- 1
thr_per_chain <- 4

# Compile non-parallel
mod <- cmdstan_model("models/wtf.stan") 

fit <- mod$sample(
  data = data_list,
  iter_warmup = 0, 
  iter_sampling = 1,
  chains = 1, 
  threads_per_chain = thr_per_chain,
  init = initfn,
  parallel_chains = 1,
  fixed_param = TRUE
)

# Works...
fit$draws("log_lik", format = "matrix")

# Compile parallel
modp <- cmdstan_model("models/wtf.stan", force_recompile = TRUE,
                      cpp_options = list(stan_threads = TRUE)) 

fitp <- modp$sample(
  data = data_list,
  iter_warmup = 0, 
  iter_sampling = 1,
  chains = 1, 
  threads_per_chain = thr_per_chain,
  init = initfn,
  parallel_chains = 1,
  fixed_param = TRUE
)

# Doesn't work...
fitp$draws("log_lik", format = "matrix")
