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
  // vector[61] D;
  vector[61] E;
  matrix[61, 61] Va;
  matrix[61, 61] m_AB;
  matrix[61, 61] A;
  real m_AA;
  real sqp;
  real meanrate = 0;
  // real mABs;
  real scale;
  vector[l] likpos = rep_vector(0, l);
  real lik = 0; // log(1)
  real likposanc;
  real cti;
  real phi;
  real N;
  real muti;
  real ttheta;
  // int pos = 1;
  
  
  A = build_A(1, kappa, 1, pi_eq);
  
  for(i in 1:61){
    meanrate -= pi_eq[i] * A[i, i];
  }
  scale = (theta / 2.0) / meanrate;
  // Calculate substitution rate matrix
  // mutmat = PDRM(mu, kap, om, pi_eq);
  mutmat = build_A(1, kappa, omega, pi_eq);
  
  // Eigenvectors/values of substitution rate matrix
  V = eigenvectors_sym(mutmat);
  // D = 1 / (1 - eigenvalues_sym(mutmat));
  E = 1 / (1 - 2 * scale * eigenvalues_sym(mutmat));
  // lambda/(lambda-2.0*scale*SymEigenval[i])
  V_inv = diag_post_multiply(V, E);
  
  // Create m_AB for each ancestral codon
  for(i in 1:61) {
    Va = rep_matrix(row(V, i), 61);
    m_AB[i, ] = to_row_vector(rows_dot_product(Va, V_inv));
  }
  
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
  
  // Likelihood calculation
  for(pos in 1:l) {
    
    phi = 0;
    N = 0;
    
    for(i in 1:61){
      cti = X[pos, i];
      phi -= lgamma(cti + 1);
      N += cti;
    }
    
    phi += lgamma(N + 1);
    
    for(anc in 1:61){
      
      likposanc = log(pi_eq[anc]);
      ttheta = 0;
      
      for(i in 1:61){
        cti = X[pos, i];
        muti = m_AB[anc, i];
        ttheta += muti;
        
        if(i == anc){
          likposanc += lgamma(cti + muti + 1) - lgamma(muti + 1);
        }else{
          likposanc += lgamma(cti + muti) - lgamma(muti);
        }
      }
      
      likposanc += lgamma(ttheta) - lgamma(N + ttheta) - log(N + ttheta) + log(ttheta);
      likpos[pos] = sillyplus(likpos[pos], likposanc);
      // likpos[pos] += likposanc;
    }
   // likpos[pos] = sillymult(likpos[pos], phi);
   likpos[pos] += phi;
   lik += likpos[pos];
  }
}


