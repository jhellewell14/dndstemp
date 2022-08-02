#include functions.stan

data {
  int <lower = 1> l; // gene length
  int <lower = 1> n; // number of samples
  array[l, 61] int X; // codon frequencies at each locus in the gene
  row_vector[61] pi_eq; // equilibrium codon frequencies
  array[l] int N;
  int grainsize;
} 

parameters {
  real om;
  real kap;
  real th;
}

transformed parameters {
  real m_AA;
  matrix[61, 61] m_AB;
  matrix[61, 61] mutmat;
  matrix[61, 61] A;
  matrix[61, 61] V;
  matrix[61, 61] Va;
  vector[61] D;
  vector[61] E;
  matrix[61, 61] V_inv;
  real sqp;
  real meanrate = 0;
  real scale;
  real omega = exp(om);
  real kappa = exp(kap);
  real theta = exp(th);
  
  profile("mats"){
  // Build mutation rate matrix under neutrality
  // https://github.com/danny-wilson/genomegaMap/blob/659802ef273af54e1b69af77dcfae210eb250e3f/src/genomegaMap/Utilities/mutation.cpp#L1448
  A = build_A(1, kappa, 1, pi_eq);
  
  // Calculate mean rate
  // https://github.com/danny-wilson/genomegaMap/blob/659802ef273af54e1b69af77dcfae210eb250e3f/src/genomegaMap/Utilities/mutation.cpp#L1451
  for(i in 1:61){
    meanrate -= pi_eq[i] * A[i, i];
  }
  scale = (theta / 2.0) / meanrate;
  
  // Calculate substitution rate matrix again
  // not under neutrality
  // https://github.com/danny-wilson/genomegaMap/blob/659802ef273af54e1b69af77dcfae210eb250e3f/src/genomegaMap/Utilities/mutation.cpp#L1458
  // mutmat = PDRM(mu, kap, om, pi_eq); (old)
  mutmat = build_A(1, kappa, omega, pi_eq);
  }

  // Eigenvectors/values of substitution rate matrix
  V = eigenvectors_sym(mutmat);
  // D = 1 / (1 - eigenvalues_sym(mutmat)); (old)
  
  // These few lines correspond to
  // https://github.com/danny-wilson/genomegaMap/blob/659802ef273af54e1b69af77dcfae210eb250e3f/src/genomegaMap/Utilities/mutation.cpp#L221-L232
  E = 1 / (1 - 2 * scale * eigenvalues_sym(mutmat));
  V_inv = diag_post_multiply(V, E);
  // Create m_AB for each ancestral codon
  for(i in 1:61) {
    Va = rep_matrix(row(V, i), 61);
    m_AB[i, ] = to_row_vector(rows_dot_product(Va, V_inv));
  }
  
  // Add equilibrium frequencies
  // https://github.com/danny-wilson/genomegaMap/blob/659802ef273af54e1b69af77dcfae210eb250e3f/src/genomegaMap/Utilities/mutation.cpp#L233
  for(i in 1:61){
    sqp = sqrt(pi_eq[i]);
    m_AB[i, ] /= sqp;
    m_AB[, i] *= sqp;
  }
  
  // Normalise by m_AA
  // https://github.com/danny-wilson/genomegaMap/blob/659802ef273af54e1b69af77dcfae210eb250e3f/src/genomegaMap/Transformations/NY98_PDRM.cpp#L97-L109
  for(i in 1:61){
    m_AA = m_AB[i, i];
    for(j in 1:61){
      if(j != i){
        m_AB[i, j] /= m_AA;
      }
      if(m_AB[i, j] < 0){
        m_AB[i, j] = 1.0e-06;
      }
    }
    m_AB[i, i] = 1.0e-06;
  }
  
}

model { 

  target += reduce_sum(partial_sum_constant, N,
  grainsize,
  m_AB,
  X);

  kap ~ normal(0, 1);
  th ~ normal(0, 1);
  om ~ normal(0, 1);
}
