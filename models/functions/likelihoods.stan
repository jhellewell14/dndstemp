

// vector f(vector beta, vector gam,
// data array[] real x_r, data array[] int x_i){
//   int n = size (x_i);
//   int l = n / 61;
//   array[l] matrix[61, 61] obs_mat;
//   array[l] vector[61] obs_vec;
//   vector[61] pi_eq = to_vector(x_r);
//   matrix[61, 61] pimat = diag_matrix(sqrt(to_vector(pi_eq)));
//   matrix[61, 61] pimatinv = diag_matrix(inv(sqrt(to_vector(pi_eq))));
//   vector[61] ones = rep_vector(1, 61);
//   matrix[61, 61] pimult;
//   matrix[61, 61] mutmat;
//   matrix[61, 61] V;
//   matrix[61, 61] V_inv;
//   vector[61] E;
//   matrix[61, 61] Va;
//   matrix[61, 61] m_AB;
//   matrix[61, 61] A;
//   real m_AA;
//   real sqp;
//   real meanrate = 0;
//   real scale;
//   real lik = 0; // log(1)
//   vector[61] likposanc;
//   real phi;
//   int Np;
//   matrix[61, 61] muti;
//   matrix[61, 61] lgmuti;
//   matrix[61, 61] gam_mat;
//   vector[61] ttheta;
//   vector[61] ltheta;
//   vector[61] lgtheta;
//   matrix[61, 61] mti;
//   vector[61] x;
//   vector[61] poslp;
//   array[l] int N;
//   vector[l] lv;
//   vector[61] lp = to_vector(log(pi_eq));
//   vector[l] out;
// 
//   // Reconstruct data
//   real theta = beta[1];
//   real kappa = beta[2];
//   real omega = beta[3];
// 
//   for(j in 1:61){
//     for(i in 1:61){
//       pimult[i, j] = sqrt(pi_eq[j] / pi_eq[i]);
//     }
//   }
// 
//   for(i in 1:l){
//     int j = 1 + ((i - 1) * 61);
//     int k = j + 60;
//     obs_vec[i] = to_vector(x_i[j:k]);
//     obs_mat[i] = rep_matrix(to_row_vector(obs_vec[i]), 61);
//     N[i] = sum(x_i[j:k]);
//   }
// 
// 
//   // Find mean mutation rate under neutrality
//   A = build_A(kappa, 1, pimat, pimult);
//   meanrate = - dot_product(pi_eq, diagonal(A));
//   scale = (theta / 2.0) / meanrate;
// 
//   // Calculate substitution rate matrix not under neutrality
//   // here you can multiply the relevant elements by new omega value and recalculate the diagonal
//   mutmat = update_A(A, omega, pimult);
// 
//   // Eigenvectors/values of substitution rate matrix
//   V = eigenvectors_sym(mutmat);
//   E = 1 / (1 - 2 * scale * eigenvalues_sym(mutmat));
//   V_inv = diag_post_multiply(V, E);
// 
//   // Create m_AB for each ancestral codon
//   for(i in 1:61) {
//     Va = rep_matrix(row(V, i), 61);
//     m_AB[, i] = rows_dot_product(Va, V_inv);
//   }
// 
//   // Multiply by equilibrium frequencies
//   // m_AB[, i] /= sqrt(pi_eq[i]);
//   // m_AB[i, ] *= sqrt(pi_eq[i]);
//   m_AB = (m_AB' * pimatinv)' * pimat;
// 
//   // Normalise - m_AB / m_AA
//   for(i in 1:61){
//     // sqp = sqrt(pi_eq[i]);
//     // m_AB[, i] /= sqp;
//     // m_AB[i, ] *= sqp;
//     m_AB[, i] /= m_AB[i, i];
//     m_AB[i, i] = 1.0e-06; // This happens in the C code but not mentioned elsewhere
//     for(j in 1:61){if(m_AB[i, j] < 0) m_AB[i, j] = 1.0e-06;}
//   }
// 
//   // Writing to columns was faster so now we transpose
//   m_AB = m_AB';
// 
//   // Likelihood calculation
// 
//   // Parts shared over all positions (at least while omega is fixed)
//   muti = add_diag(m_AB, 1);
//   lgmuti = lgamma(muti);
//   ttheta = m_AB * ones;
//   ltheta = log(ttheta);
//   lgtheta = lgamma(ttheta);
// 
//   for(pos in 1:l) {
// 
//     // Calculate parts same for all ancestors
//     Np = N[pos];
//     phi = lgamma(Np + 1) - sum(lgamma(obs_vec[pos] + 1));
//     poslp = lgtheta - lgamma(Np + ttheta) - log(Np + ttheta) + ltheta;
// 
//     gam_mat = lgamma(obs_mat[pos] + muti) - lgmuti;
// 
//     likposanc = lp;
//     likposanc += gam_mat * ones;
//     likposanc += poslp + phi;
//     lv[pos] = log_sum_exp(likposanc);
//     lik += log_sum_exp(likposanc);
//   }
// 
//   out = lv;
//   // out[1] = lik;
//   return out;
// }

real rs(array[] matrix X_slice,
int start, int end,
matrix A,
matrix pimult,
matrix pimatinv,
matrix pimat,
real scale,
array[] int N,
vector omega,
vector lp,
vector phi) {
  
  matrix[61, 61] mutmat, V, V_inv, m_AB, Va, obs_mat, muti, lgmuti, gam_mat;
  vector[61] E, obs_vec, ttheta, ltheta, lgtheta, poslp, likposanc;
  vector[61] ones = rep_vector(1, 61);
  real lik = 0;
  real Np;
  int z = 1;
  
  for(pos in start:end){
    obs_vec = to_vector(row(X_slice[z], 1));
    obs_mat = X_slice[z];
    z += 1;
    
    // Calculate substitution rate matrix not under neutrality
    // here you can multiply the relevant elements by new omega value and recalculate the diagonal
    mutmat = update_A(A, omega[pos], pimult);
    
    // Eigenvectors/values of substitution rate matrix
    V = eigenvectors_sym(mutmat);
    E = 1 / (1 - 2 * scale * eigenvalues_sym(mutmat));
    V_inv = diag_post_multiply(V, E);
    
    // Create m_AB for each ancestral codon
    for(i in 1:61) {
      Va = rep_matrix(row(V, i), 61);
      m_AB[, i] = rows_dot_product(Va, V_inv);
    }
    
    // Multiply by equilibrium frequencies
    m_AB = (m_AB' * pimatinv)' * pimat;
    
    // Normalise - m_AB / m_AA
    for(i in 1:61){
      m_AB[, i] /= m_AB[i, i];
      m_AB[i, i] = 1.0e-06; // This happens in the C code but not mentioned elsewhere
      for(j in 1:61){if(m_AB[i, j] < 0) m_AB[i, j] = 1.0e-06;}
    }
    
    // Writing to columns was faster so now we transpose
    m_AB = m_AB';
    
    // Likelihood calculation
    muti = add_diag(m_AB, 1);
    lgmuti = lgamma(muti);
    ttheta = m_AB * ones;
    ltheta = log(ttheta);
    lgtheta = lgamma(ttheta);
    
    // Calculate parts same for all ancestors
    Np = N[pos];
    // phi = lgamma(Np + 1) - sum(lgamma(obs_vec + 1));
    poslp = lgtheta - lgamma(Np + ttheta) - log(Np + ttheta) + ltheta;
    
    // Calculates parts different for each ancestor
    gam_mat = lgamma(obs_mat + muti) - lgmuti;
    
    // Sum together
    likposanc = lp;
    likposanc += gam_mat * ones;
    likposanc += poslp + phi[pos];
    lik += log_sum_exp(likposanc);
  }
  return lik;
}


real rs_GTR(array[] matrix X_slice,
int start, int end,
matrix A,
matrix pimult,
matrix pimatinv,
matrix pimat,
real scale,
array[] int N,
vector omega,
real alpha,
real beta,
real gamma,
real delta,
real epsilon,
real eta,
vector lp,
vector phi) {
  
  matrix[61, 61] mutmat, V, V_inv, m_AB, Va, obs_mat, muti, lgmuti, gam_mat;
  vector[61] E, obs_vec, ttheta, ltheta, lgtheta, poslp, likposanc;
  vector[61] ones = rep_vector(1, 61);
  real lik = 0;
  real Np;
  int z = 1;
  
  for(pos in start:end){
    obs_vec = to_vector(row(X_slice[z], 1));
    obs_mat = X_slice[z];
    z += 1;
    
    // Calculate substitution rate matrix not under neutrality
    // here you can multiply the relevant elements by new omega value and recalculate the diagonal
    mutmat = build_GTR(alpha, beta, gamma, delta, epsilon, eta, omega[pos], pimat, pimult);
    
    // Eigenvectors/values of substitution rate matrix
    V = eigenvectors_sym(mutmat);
    E = 1 / (1 - 2 * scale * eigenvalues_sym(mutmat));
    V_inv = diag_post_multiply(V, E);
    
    // Create m_AB for each ancestral codon
    for(i in 1:61) {
      Va = rep_matrix(row(V, i), 61);
      m_AB[, i] = rows_dot_product(Va, V_inv);
    }
    
    // Multiply by equilibrium frequencies
    m_AB = (m_AB' * pimatinv)' * pimat;
    
    // Normalise - m_AB / m_AA
    for(i in 1:61){
      m_AB[, i] /= m_AB[i, i];
      m_AB[i, i] = 1.0e-06; // This happens in the C code but not mentioned elsewhere
      for(j in 1:61){if(m_AB[i, j] < 0) m_AB[i, j] = 1.0e-06;}
    }
    
    // Writing to columns was faster so now we transpose
    m_AB = m_AB';
    
    // Likelihood calculation
    muti = add_diag(m_AB, 1);
    lgmuti = lgamma(muti);
    ttheta = m_AB * ones;
    ltheta = log(ttheta);
    lgtheta = lgamma(ttheta);
    
    // Calculate parts same for all ancestors
    Np = N[pos];
    // phi = lgamma(Np + 1) - sum(lgamma(obs_vec + 1));
    poslp = lgtheta - lgamma(Np + ttheta) - log(Np + ttheta) + ltheta;
    
    // Calculates parts different for each ancestor
    gam_mat = lgamma(obs_mat + muti) - lgmuti;
    
    // Sum together
    likposanc = lp;
    likposanc += gam_mat * ones;
    likposanc += poslp + phi[pos];
    lik += log_sum_exp(likposanc);
  }
  return lik;
}


