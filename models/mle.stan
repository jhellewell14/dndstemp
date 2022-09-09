functions{
#include functions/PDRM_mat.stan
#include functions/deprecated.stan
}

data {
  row_vector[61] pi_eq;
  int <lower = 1> l; // gene length
  array[l, 61] int X; // codon frequencies at each locus in the gene
}

transformed data {
  matrix[61, 61] pimat = diag_matrix(sqrt(to_vector(pi_eq)));
  matrix[61, 61] pimatinv = diag_matrix(inv(sqrt(to_vector(pi_eq))));
  matrix[61, 61] pimult;
  
  for(j in 1:61){
    for(i in 1:61){
      pimult[i, j] = sqrt(pi_eq[j] / pi_eq[i]);
    }
  }
}
 
parameters {
  real <lower = 0> kappa;
  real <lower = 0> omega;
  real <lower = 0> theta;
}

model {
  matrix[61, 61] mutmat;
  matrix[61, 61] V;
  matrix[61, 61] V_inv;
  vector[61] E;
  vector[61] Ve;
  matrix[61, 61] Va;
  matrix[61, 61] m_AB;
  matrix[61, 61] A;
  real m_AA;
  real sqp;
  real meanrate = 0;
  real scale;
  vector[l] likpos = rep_vector(0, l);
  real likposanc;
  real cti;
  real phi;
  real N;
  real muti;
  real ttheta;
  
  A = build_A(kappa, 1, pimat, pimult);
  
  for(i in 1:61){
    meanrate -= pi_eq[i] * A[i, i];
  }
  scale = (theta / 2.0) / meanrate;
  // Calculate substitution rate matrix
  // mutmat = PDRM(mu, kap, om, pi_eq);
  mutmat = update_A(A, omega, pimult); 
  
  // Eigenvectors/values of substitution rate matrix
  V = eigenvectors_sym(mutmat);
  Ve = eigenvalues_sym(mutmat);
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
    }

   likpos[pos] += phi;
   // lik += likpos[pos];
   target += likpos[pos];
  }
}
