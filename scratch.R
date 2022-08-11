library(posterior)
library(bayesplot)
library(data.table)
library(magrittr)
library(ggplot2)
library(cmdstanr)
library(patchwork)
library(truncnorm)
library(seqinr)
source("generate_data.R")

# Read alignment
data <- seqinr::read.alignment(file = "fits/porB/porB3.carriage.noindels.txt", format = "fasta")
# Create list of data for stan
# data_list$X is the codon counts at each location
data_list <- generate_data(data)
data_list$X <- data_list$X[1:100, ]
data_list$N <- data_list$N[1:100]
data_list$l <- 100
data_list$om_sd_prior <- 6
data_list$theta_prior <- 2

# Stan parameters
nchains <- 1
cores <- (parallel::detectCores() - 2)
thr_per_chain <- floor(cores / nchains)
Sys.setenv(STAN_NUM_THREADS = 8)

# Compile model
mod <- cmdstan_model("models/o_hier_mk_fixed.stan", 
                     cpp_options = list(stan_threads = TRUE))

# Fit model
fit <- mod$sample(
  data = data_list,
  iter_warmup = 500, 
  iter_sampling = 250,
  chains = 2, 
  threads_per_chain = thr_per_chain,
  init = initfn,
  parallel_chains = 1,
  refresh = 50
)

# Output
fit$summary(c("theta[1]", "kappa[1]"))
fit$summary()
res <- fit$draws(format = "array")

mcmc_intervals(res, regex_pars = "omega\\[") +
  geom_vline(xintercept = median(as.vector(fit$draws("omega", format = "matrix"))), col = "red", lty = 2)
  # geom_vline(xintercept = exp(median(as.vector(fit$draws("om_mean", format = "matrix")))))
mcmc_hist(res, pars = c("th", "kap"))
lp_ncp <- log_posterior(fit)
np_ncp <- nuts_params(fit)


mcmc_trace(res, pars = "th", np = np_ncp)
# mcmc_parcoord(res, pars = c("om_mean", "om_sd",
#                                   "kap", "th", "om_raw[1]"), np = np_ncp)

mcmc_scatter(res, pars = c("kap", "th"), np = np_ncp,
             , transformations = list())
mcmc_pairs(res, pars = c("kap", "th", "om_raw[1]",
                         "om_raw[9]", "om_raw[15]",
                         "om_raw[12]", "om_raw[5]"), lp = lp_ncp,
           np = np_ncp, transformations = list("kap" = exp, "th" = exp))


##### IGNORE AFTER HERE UNTIL MODEL WORKS




DT <- data.table::rbindlist(list(extract_res(fit = fit, mod_name = "fix")))
DT
hist(exp(as.vector(fit$draws("kap"))))
# DT <- data.table::rbindlist(lapply(list.files(path = "fits/penA/",pattern = "*.RDS"), function(x){
#   extract_res(fit = readRDS(paste0("fits/penA/", x)), 
#               mod_name = strsplit(x, split = "[.]")[[1]][1])
# }))

fit <- readRDS("fits/porB/o_hier_mk_fixed.RDS")


line_df <- data.frame(par = factor("omega", levels = c("omega", "kappa", "mu")), 
                      z = 1)

DT[] %>%
  ggplot(aes(x = siteno, y = med, ymin = lq, ymax = uq, col = mod)) +
  geom_point(position = position_dodge(width = 1)) +
  geom_errorbar(position = position_dodge(width = 1)) +
  facet_wrap(par ~ ., ncol = 1, scales = "free_y") +
  scale_y_log10() + 
  geom_hline(data = line_df, aes(yintercept = z), lty = 2)

DT[par == "kappa"]

plot_dnds(data_list$X)
