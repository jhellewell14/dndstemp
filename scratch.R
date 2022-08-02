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

# Stan parameters
nchains <- 1
cores <- (parallel::detectCores() - 2)
thr_per_chain <- floor(cores / nchains)
Sys.setenv(STAN_NUM_THREADS = 8)

# Compile model
mod <- cmdstan_model("models/constant_o.stan", 
                     cpp_options = list(stan_threads = TRUE))

# Fit model
fit <- mod$sample(
  data = data_list,
  iter_warmup = 500, 
  iter_sampling = 500,
  chains = nchains, 
  threads_per_chain = thr_per_chain,
  init = initfn,
  parallel_chains = 4,
  refresh = 100 # print update every 100 iters
)

# Output
fit$summary(c("theta", "omega", "kappa"))


##### IGNORE AFTER HERE UNTIL MODEL WORKS


DT <- data.table::rbindlist(list(extract_res(fit = fit, mod_name = "fix")))
DT
hist(as.vector(fit$draws("kappa")))
# DT <- data.table::rbindlist(lapply(list.files(path = "fits/penA/",pattern = "*.RDS"), function(x){
#   extract_res(fit = readRDS(paste0("fits/penA/", x)), 
#               mod_name = strsplit(x, split = "[.]")[[1]][1])
# }))

# fit <- readRDS("fits/porB/o_hier_mk_fixed.RDS")


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
