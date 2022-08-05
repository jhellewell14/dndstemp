#include functions.stan

data {
  row_vector[61] pi_eq;
  int <lower = 1> l; // gene length
  // int loc;
  array[l, 61] int X; // codon frequencies at each locus in the gene
  real <lower = 0> mu;
  real <lower = 0> theta;
  // real <lower = 0> theta;
  real <lower = 0> kappa;
  real <lower = 0> omega;
  array [l] int N;
}

transformed data{
  vector[61] lp = to_vector(log(pi_eq));
  array[l] matrix[61, 61] obs_mat;
  array[l] vector[61] obs_vec;
  vector[61] ones = rep_vector(1, 61);
  
  for(i in 1:l){
    obs_mat[i] = rep_matrix(to_row_vector(X[i, ]), 61);
    obs_vec[i] = to_vector(X[i, ]);
  }
  
}

parameters {
  // real <lower = 0> mu;
  // real <lower = 0> kappa;
  // real <lower = 0> omega;
  // real <lower = 0> theta;
}

model {
  
}

generated quantities {
  matrix[61, 61] mutmat = PDRM(mu, kappa, omega, pi_eq);
  matrix[61, 61] V;
  matrix[61, 61] V_inv;
  vector[61] E;
  matrix[61, 61] Va;
  matrix[61, 61] m_AB;
  matrix[61, 61] A;
  real m_AA;
  real sqp;
  real meanrate = 0;
  real scale;
  real lik = 0; // log(1)
  vector[61] likposanc;
  real phi;
  int Np;
  matrix[61, 61] muti;
  matrix[61, 61] lgmuti;
  matrix[61, 61] gam_mat;
  vector[61] ttheta;
  vector[61] ltheta;
  vector[61] lgtheta;
  matrix[61, 61] mti;
  vector[61] x;
  vector[61] poslp;
  
  // Find mean mutation rate under neutrality
  A = build_A(1, kappa, 1, pi_eq);
  meanrate = - dot_product(pi_eq, diagonal(A));
  scale = (theta / 2.0) / meanrate;
  
  // Calculate substitution rate matrix not under neutrality
  mutmat = build_A(1, kappa, omega, pi_eq); 
  
  // Eigenvectors/values of substitution rate matrix
  V = eigenvectors_sym(mutmat);
  E = 1 / (1 - 2 * scale * eigenvalues_sym(mutmat));
  V_inv = diag_post_multiply(V, E);
  
  // Create m_AB for each ancestral codon
  for(i in 1:61) {
    Va = rep_matrix(row(V, i), 61);
    m_AB[, i] = rows_dot_product(Va, V_inv);
  }
  
  m_AB = m_AB';
  
  // Add equilibrium frequencies
  for(i in 1:61){
    sqp = sqrt(pi_eq[i]);
    m_AB[i, ] /= sqp;
    m_AB[, i] *= sqp;
  }
  
  // Normalise by m_AA
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
  
  muti = add_diag(m_AB, 1);
  lgmuti = lgamma(muti);
  ttheta = m_AB * ones;
  ltheta = log(ttheta);
  lgtheta = lgamma(ttheta);
  // Likelihood calculation
  for(pos in 1:l) {
    
    // Calculate parts same for all ancestors
    Np = N[pos];
    phi = lgamma(Np + 1) - sum(lgamma(obs_vec[pos] + 1));
    
    poslp = lgtheta - lgamma(Np + ttheta) - log(Np + ttheta) + ltheta;
    
    gam_mat = lgamma(obs_mat[pos] + muti) - lgmuti;
    
    likposanc = lp;
    likposanc += gam_mat * ones;
    likposanc += poslp + phi;

    lik += log_sum_exp(likposanc);
  }
}


