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
# data <- seqinr::read.alignment(file = "~/Downloads/for_Joel/SPARC1_CLS00381.dna.aln", format = "fasta")
# Create list of data for stan
# data_list$X is the codon counts at each location
data_list <- generate_data(data)

# Stan parameters
nchains <- 1
cores <- (parallel::detectCores() - 2)
thr_per_chain <- floor(cores / nchains)
Sys.setenv(STAN_NUM_THREADS = 8)

# mle <- cmdstan_model("models/mle.stan")

# mle_fit <- mle$optimize(data = data_list,
#                         init = list(list("omega" = 0.5, "theta" = 0.5, "kappa" = 1,
#                                          "alpha" = 0.5, "beta" = 0.5, "gamma" = 0.5,
#                                          "delta" = 0.5, "epsilon" = 0.5, "eta" = 0.5)))
# 
# mle_fit$summary()

# data_list$kp_prior <- as.numeric(mle_fit$draws("kappa"))
# data_list$th_prior <- as.numeric(mle_fit$draws("theta"))
# data_list$om_prior <- as.numeric(mle_fit$draws("omega"))
# data_list$om_sd_prior <- 1
# k <- 100
# data_list$X <- data_list$X[1:k, ]
# data_list$N <- data_list$N[1:k]
# data_list$l <- k

# Compile model
mod <- cmdstan_model("models/mv.stan", 
                     cpp_options = list(stan_threads = TRUE), 
                     stanc_options = list("O1"))

# Fit model
fit <- mod$sample(
  data = data_list,
  iter_warmup = 500, 
  iter_sampling = 250,
  chains = 1, 
  threads_per_chain = thr_per_chain,
  init = initfn,
  parallel_chains = 1,
  refresh = 50
)

ff <- mod$variational(data = data_list, init = initfn, threads = 8)



# Output
fit$summary(c("theta[1]", "kappa[1]"))
res <- ff$draws(format = "array")


mcmc_intervals(res, regex_pars = "omega\\[") +
  geom_vline(xintercept = exp(median(as.vector(ff$draws("om_mean", format = "matrix")))), col = "red", lty = 2)
  
lp_ncp <- log_posterior(fit)
np_ncp <- nuts_params(fit)

mcmc_scatter(res, pars = c("om_raw[97]", "lp__"), np = np_ncp)

mcmc_pairs(res, pars = c("kap", "th", "om_mean", paste0("om_raw[",96:100,"]")), lp = lp_ncp,
           np = np_ncp)


##### IGNORE AFTER HERE UNTIL MODEL WORKS




DT <- data.table::rbindlist(list(extract_res(fit = fit, mod_name = "fix")))
DT
hist(exp(as.vector(fit$draws("kap"))))
# DT <- data.table::rbindlist(lapply(list.files(path = "fits/penA/",pattern = "*.RDS"), function(x){
#   extract_res(fit = readRDS(paste0("fits/penA/", x)), 
#               mod_name = strsplit(x, split = "[.]")[[1]][1])
# }))

fit <- readRDS("fits/penA/all_idd.RDS")
fit <- readRDS("fits/pneumo/SPARC1_CLS01692/hier_exp2.RDS")
fit$diagnostic_summary()
fitgtr$diagnostic_summary()



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
