vector f2(vector beta_vec, vector gam,
data array[] real x_r, data array[] int x_i){
  int n = size (x_i);
  int shard_ends = x_i[n];
  int shard_starts = x_i[n - 1];
  int K = x_i[n - 2];
  int l = x_i[n - 3];
  array[l] matrix[61, 61] obs_mat;
  vector[61] obs_vec;
  vector[61] lp = to_vector(x_r[1:61]);
  vector[61] ones = rep_vector(1, 61);
  vector[61] likposanc;
  vector[l] phi = to_vector(x_r[(61 + shard_starts):(61 + shard_ends)]);
  int Np;
  matrix[61, 61] muti;
  matrix[61, 61] m_AB = to_matrix(beta_vec, 61, 61, 0); // this fils in matrix in row order
  matrix[61, 61] lgmuti;
  matrix[61, 61] gam_mat;
  vector[61] ttheta;
  vector[61] ltheta;
  vector[61] lgtheta;
  vector[61] poslp;
  array[l] int N;
  vector[l] lv;
  vector[l] out;
  
  for(i in 1:l){
    int j = 1 + ((i - 1) * 61);
    int k = j + 60;
    obs_vec = to_vector(x_i[j:k]);
    obs_mat[i] = rep_matrix(to_row_vector(obs_vec), 61);
    N[i] = sum(x_i[j:k]);
    // m_AB[i, ] = to_row_vector(beta_vec[((i - 1) * 61 + 1):((i - 1) * 61 + 61)]);
  }
  
  // // Likelihood calculation
  // // Parts shared over all positions (at least while omega is fixed)
  muti = add_diag(m_AB, 1);
  lgmuti = lgamma(muti);
  ttheta = m_AB * ones;
  ltheta = log(ttheta);
  lgtheta = lgamma(ttheta);

  for(pos in 1:l) {

    // Calculate parts same for all ancestors
    Np = N[pos];
    poslp = lgtheta - lgamma(Np + ttheta) - log(Np + ttheta) + ltheta;

    gam_mat = lgamma(obs_mat[pos] + muti) - lgmuti;

    likposanc = lp;
    likposanc += gam_mat * ones;
    likposanc += poslp + phi[pos];
    // Indices makes sure that for each location the likelihoods for each mixture
    // component are next to each other
    print(likposanc);
    lv[pos] = log_sum_exp(likposanc);
  }


  out = lv;
  return out;
}