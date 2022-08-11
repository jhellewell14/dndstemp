functions {
#include functions/GTR_mat.stan
#include functions/PDRM_mat.stan
#include functions/likelihoods.stan
}

data {
  row_vector[61] pi_eq;
  int <lower = 1> l; // gene length
  array[l, 61] int X; // codon frequencies at each locus in the gene
  array [l] int N;
  real om_sd_prior;
  real theta_prior;
} 

transformed data{
  vector[61] lp = to_vector(log(pi_eq));
  array[l] matrix[61, 61] obs_mat;
  matrix[61, 61] pimat = diag_matrix(sqrt(to_vector(pi_eq)));
  matrix[61, 61] pimatinv = diag_matrix(inv(sqrt(to_vector(pi_eq))));
  matrix[61, 61] pimult;
  int grainsize = 1;
  
  for(j in 1:61){
    for(i in 1:61){
      pimult[i, j] = sqrt(pi_eq[j] / pi_eq[i]);
    }
  }

  for(i in 1:l){
    obs_mat[i] = rep_matrix(to_row_vector(X[i, ]), 61);
  }
  
}

parameters {
  real kap;
  // real om_mean;
  // real <lower = 0> om_sd;
  vector [l] om_raw;
  real th;
} 

transformed parameters {
  // Find mean mutation rate under neutrality 
  // (same for all locations because kappa shared across locations)
  vector[l] kappa = rep_vector(exp(kap), l);
  vector[l] theta = rep_vector(exp(th), l);
  // vector[l] omega = exp(om_mean + om_raw * om_sd);
  vector[l] omega = exp(om_raw);
  matrix[61, 61] A = build_A(kap, 1, pimat, pimult);
  real meanrate = - dot_product(pi_eq, diagonal(A));
  real scale = (th / 2.0) / meanrate;
}

model {
  // Priors
  // om_mean ~ normal(0, 1);
  om_raw ~ normal(0, 1);
  // om_sd ~ exponential(om_sd_prior);
  kap ~ normal(0, 1);
  th ~ normal(0, 1);
  // th ~ exponential(theta_prior);

  // Likelihood
  target += reduce_sum(rs,
  obs_mat,
  grainsize,
  A,
  pimult,
  pimatinv,
  pimat,
  scale,
  N,
  omega,
  lp);
}
