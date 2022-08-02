functions {
  
  // Copied from Danny's code
  // https://github.com/danny-wilson/genomegaMap/blob/659802ef273af54e1b69af77dcfae210eb250e3f/src/genomegaMap/Utilities/mutation.cpp#L1143
  matrix build_A(real mu, real kappa, real omega, row_vector pi_eq) {
    matrix[61,61] M = rep_matrix(0.,61,61);
    real sqp;
    
    M[1,2] = kappa;
    M[1,3] = omega;
    M[1,4] = omega;
    M[1,5] = omega * kappa;
    M[1,9] = omega;
    M[1,11] = omega;
    M[1,14] = omega * kappa;
    M[1,30] = omega;
    M[1,46] = omega;
    M[2,1] = kappa;
    M[2,3] = omega;
    M[2,4] = omega;
    M[2,6] = omega * kappa;
    M[2,10] = omega;
    M[2,12] = omega;
    M[2,15] = omega * kappa;
    M[2,31] = omega;
    M[2,47] = omega;
    M[3,1] = omega;
    M[3,2] = omega;
    M[3,4] = kappa;
    M[3,7] = omega * kappa;
    M[3,16] = kappa;
    M[3,32] = omega;
    M[3,48] = omega;
    M[4,1] = omega;
    M[4,2] = omega;
    M[4,3] = kappa;
    M[4,8] = omega * kappa;
    M[4,13] = omega;
    M[4,17] = kappa;
    M[4,33] = omega;
    M[4,49] = omega;
    M[5,1] = omega * kappa;
    M[5,6] = kappa;
    M[5,7] = 1;
    M[5,8] = 1;
    M[5,9] = omega;
    M[5,11] = omega;
    M[5,18] = omega * kappa;
    M[5,34] = omega;
    M[5,50] = omega;
    M[6,2] = omega * kappa;
    M[6,5] = kappa;
    M[6,7] = 1;
    M[6,8] = 1;
    M[6,10] = omega;
    M[6,12] = omega;
    M[6,19] = omega * kappa;
    M[6,35] = omega;
    M[6,51] = omega;
    M[7,3] = omega * kappa;
    M[7,5] = 1;
    M[7,6] = 1;
    M[7,8] = kappa;
    M[7,20] = omega * kappa;
    M[7,36] = omega;
    M[7,52] = omega;
    M[8,4] = omega * kappa;
    M[8,5] = 1;
    M[8,6] = 1;
    M[8,7] = kappa;
    M[8,13] = omega;
    M[8,21] = omega * kappa;
    M[8,37] = omega;
    M[8,53] = omega;
    M[9,1] = omega;
    M[9,5] = omega;
    M[9,10] = kappa;
    M[9,11] = omega * kappa;
    M[9,22] = omega * kappa;
    M[9,38] = omega;
    M[9,54] = omega;
    M[10,2] = omega;
    M[10,6] = omega;
    M[10,9] = kappa;
    M[10,12] = omega * kappa;
    M[10,23] = omega * kappa;
    M[10,39] = omega;
    M[10,55] = omega;
    M[11,1] = omega;
    M[11,5] = omega;
    M[11,9] = omega * kappa;
    M[11,12] = kappa;
    M[11,13] = omega;
    M[11,26] = omega * kappa;
    M[11,42] = omega;
    M[11,58] = omega;
    M[12,2] = omega;
    M[12,6] = omega;
    M[12,10] = omega * kappa;
    M[12,11] = kappa;
    M[12,13] = omega;
    M[12,27] = omega * kappa;
    M[12,43] = omega;
    M[12,59] = omega;
    M[13,4] = omega;
    M[13,8] = omega;
    M[13,11] = omega;
    M[13,12] = omega;
    M[13,29] = omega * kappa;
    M[13,45] = omega;
    M[13,61] = omega;
    M[14,1] = omega * kappa;
    M[14,15] = kappa;
    M[14,16] = 1;
    M[14,17] = 1;
    M[14,18] = omega * kappa;
    M[14,22] = omega;
    M[14,26] = omega;
    M[14,30] = omega;
    M[14,46] = omega;
    M[15,2] = omega * kappa;
    M[15,14] = kappa;
    M[15,16] = 1;
    M[15,17] = 1;
    M[15,19] = omega * kappa;
    M[15,23] = omega;
    M[15,27] = omega;
    M[15,31] = omega;
    M[15,47] = omega;
    M[16,3] = kappa;
    M[16,14] = 1;
    M[16,15] = 1;
    M[16,17] = kappa;
    M[16,20] = omega * kappa;
    M[16,24] = omega;
    M[16,28] = omega;
    M[16,32] = omega;
    M[16,48] = omega;
    M[17,4] = kappa;
    M[17,14] = 1;
    M[17,15] = 1;
    M[17,16] = kappa;
    M[17,21] = omega * kappa;
    M[17,25] = omega;
    M[17,29] = omega;
    M[17,33] = omega;
    M[17,49] = omega;
    M[18,5] = omega * kappa;
    M[18,14] = omega * kappa;
    M[18,19] = kappa;
    M[18,20] = 1;
    M[18,21] = 1;
    M[18,22] = omega;
    M[18,26] = omega;
    M[18,34] = omega;
    M[18,50] = omega;
    M[19,6] = omega * kappa;
    M[19,15] = omega * kappa;
    M[19,18] = kappa;
    M[19,20] = 1;
    M[19,21] = 1;
    M[19,23] = omega;
    M[19,27] = omega;
    M[19,35] = omega;
    M[19,51] = omega;
    M[20,7] = omega * kappa;
    M[20,16] = omega * kappa;
    M[20,18] = 1;
    M[20,19] = 1;
    M[20,21] = kappa;
    M[20,24] = omega;
    M[20,28] = omega;
    M[20,36] = omega;
    M[20,52] = omega;
    M[21,8] = omega * kappa;
    M[21,17] = omega * kappa;
    M[21,18] = 1;
    M[21,19] = 1;
    M[21,20] = kappa;
    M[21,25] = omega;
    M[21,29] = omega;
    M[21,37] = omega;
    M[21,53] = omega;
    M[22,9] = omega * kappa;
    M[22,14] = omega;
    M[22,18] = omega;
    M[22,23] = kappa;
    M[22,24] = omega;
    M[22,25] = omega;
    M[22,26] = omega * kappa;
    M[22,38] = omega;
    M[22,54] = omega;
    M[23,10] = omega * kappa;
    M[23,15] = omega;
    M[23,19] = omega;
    M[23,22] = kappa;
    M[23,24] = omega;
    M[23,25] = omega;
    M[23,27] = omega * kappa;
    M[23,39] = omega;
    M[23,55] = omega;
    M[24,16] = omega;
    M[24,20] = omega;
    M[24,22] = omega;
    M[24,23] = omega;
    M[24,25] = kappa;
    M[24,28] = omega * kappa;
    M[24,40] = omega;
    M[24,56] = omega;
    M[25,17] = omega;
    M[25,21] = omega;
    M[25,22] = omega;
    M[25,23] = omega;
    M[25,24] = kappa;
    M[25,29] = omega * kappa;
    M[25,41] = omega;
    M[25,57] = omega;
    M[26,11] = omega * kappa;
    M[26,14] = omega;
    M[26,18] = omega;
    M[26,22] = omega * kappa;
    M[26,27] = kappa;
    M[26,28] = 1;
    M[26,29] = 1;
    M[26,42] = omega;
    M[26,58] = omega;
    M[27,12] = omega * kappa;
    M[27,15] = omega;
    M[27,19] = omega;
    M[27,23] = omega * kappa;
    M[27,26] = kappa;
    M[27,28] = 1;
    M[27,29] = 1;
    M[27,43] = omega;
    M[27,59] = omega;
    M[28,16] = omega;
    M[28,20] = omega;
    M[28,24] = omega * kappa;
    M[28,26] = 1;
    M[28,27] = 1;
    M[28,29] = kappa;
    M[28,44] = 1;
    M[28,60] = omega;
    M[29,13] = omega * kappa;
    M[29,17] = omega;
    M[29,21] = omega;
    M[29,25] = omega * kappa;
    M[29,26] = 1;
    M[29,27] = 1;
    M[29,28] = kappa;
    M[29,45] = 1;
    M[29,61] = omega;
    M[30,1] = omega;
    M[30,14] = omega;
    M[30,31] = kappa;
    M[30,32] = 1;
    M[30,33] = omega;
    M[30,34] = omega * kappa;
    M[30,38] = omega;
    M[30,42] = omega;
    M[30,46] = omega * kappa;
    M[31,2] = omega;
    M[31,15] = omega;
    M[31,30] = kappa;
    M[31,32] = 1;
    M[31,33] = omega;
    M[31,35] = omega * kappa;
    M[31,39] = omega;
    M[31,43] = omega;
    M[31,47] = omega * kappa;
    M[32,3] = omega;
    M[32,16] = omega;
    M[32,30] = 1;
    M[32,31] = 1;
    M[32,33] = omega * kappa;
    M[32,36] = omega * kappa;
    M[32,40] = omega;
    M[32,44] = omega;
    M[32,48] = omega * kappa;
    M[33,4] = omega;
    M[33,17] = omega;
    M[33,30] = omega;
    M[33,31] = omega;
    M[33,32] = omega * kappa;
    M[33,37] = omega * kappa;
    M[33,41] = omega;
    M[33,45] = omega;
    M[33,49] = omega * kappa;
    M[34,5] = omega;
    M[34,18] = omega;
    M[34,30] = omega * kappa;
    M[34,35] = kappa;
    M[34,36] = 1;
    M[34,37] = 1;
    M[34,38] = omega;
    M[34,42] = omega;
    M[34,50] = omega * kappa;
    M[35,6] = omega;
    M[35,19] = omega;
    M[35,31] = omega * kappa;
    M[35,34] = kappa;
    M[35,36] = 1;
    M[35,37] = 1;
    M[35,39] = omega;
    M[35,43] = omega;
    M[35,51] = omega * kappa;
    M[36,7] = omega;
    M[36,20] = omega;
    M[36,32] = omega * kappa;
    M[36,34] = 1;
    M[36,35] = 1;
    M[36,37] = kappa;
    M[36,40] = omega;
    M[36,44] = omega;
    M[36,52] = omega * kappa;
    M[37,8] = omega;
    M[37,21] = omega;
    M[37,33] = omega * kappa;
    M[37,34] = 1;
    M[37,35] = 1;
    M[37,36] = kappa;
    M[37,41] = omega;
    M[37,45] = omega;
    M[37,53] = omega * kappa;
    M[38,9] = omega;
    M[38,22] = omega;
    M[38,30] = omega;
    M[38,34] = omega;
    M[38,39] = kappa;
    M[38,40] = omega;
    M[38,41] = omega;
    M[38,42] = omega * kappa;
    M[38,54] = omega * kappa;
    M[39,10] = omega;
    M[39,23] = omega;
    M[39,31] = omega;
    M[39,35] = omega;
    M[39,38] = kappa;
    M[39,40] = omega;
    M[39,41] = omega;
    M[39,43] = omega * kappa;
    M[39,55] = omega * kappa;
    M[40,24] = omega;
    M[40,32] = omega;
    M[40,36] = omega;
    M[40,38] = omega;
    M[40,39] = omega;
    M[40,41] = kappa;
    M[40,44] = omega * kappa;
    M[40,56] = omega * kappa;
    M[41,25] = omega;
    M[41,33] = omega;
    M[41,37] = omega;
    M[41,38] = omega;
    M[41,39] = omega;
    M[41,40] = kappa;
    M[41,45] = omega * kappa;
    M[41,57] = omega * kappa;
    M[42,11] = omega;
    M[42,26] = omega;
    M[42,30] = omega;
    M[42,34] = omega;
    M[42,38] = omega * kappa;
    M[42,43] = kappa;
    M[42,44] = omega;
    M[42,45] = omega;
    M[42,58] = omega * kappa;
    M[43,12] = omega;
    M[43,27] = omega;
    M[43,31] = omega;
    M[43,35] = omega;
    M[43,39] = omega * kappa;
    M[43,42] = kappa;
    M[43,44] = omega;
    M[43,45] = omega;
    M[43,59] = omega * kappa;
    M[44,28] = 1;
    M[44,32] = omega;
    M[44,36] = omega;
    M[44,40] = omega * kappa;
    M[44,42] = omega;
    M[44,43] = omega;
    M[44,45] = kappa;
    M[44,60] = omega * kappa;
    M[45,13] = omega;
    M[45,29] = 1;
    M[45,33] = omega;
    M[45,37] = omega;
    M[45,41] = omega * kappa;
    M[45,42] = omega;
    M[45,43] = omega;
    M[45,44] = kappa;
    M[45,61] = omega * kappa;
    M[46,1] = omega;
    M[46,14] = omega;
    M[46,30] = omega * kappa;
    M[46,47] = kappa;
    M[46,48] = 1;
    M[46,49] = 1;
    M[46,50] = omega * kappa;
    M[46,54] = omega;
    M[46,58] = omega;
    M[47,2] = omega;
    M[47,15] = omega;
    M[47,31] = omega * kappa;
    M[47,46] = kappa;
    M[47,48] = 1;
    M[47,49] = 1;
    M[47,51] = omega * kappa;
    M[47,55] = omega;
    M[47,59] = omega;
    M[48,3] = omega;
    M[48,16] = omega;
    M[48,32] = omega * kappa;
    M[48,46] = 1;
    M[48,47] = 1;
    M[48,49] = kappa;
    M[48,52] = omega * kappa;
    M[48,56] = omega;
    M[48,60] = omega;
    M[49,4] = omega;
    M[49,17] = omega;
    M[49,33] = omega * kappa;
    M[49,46] = 1;
    M[49,47] = 1;
    M[49,48] = kappa;
    M[49,53] = omega * kappa;
    M[49,57] = omega;
    M[49,61] = omega;
    M[50,5] = omega;
    M[50,18] = omega;
    M[50,34] = omega * kappa;
    M[50,46] = omega * kappa;
    M[50,51] = kappa;
    M[50,52] = 1;
    M[50,53] = 1;
    M[50,54] = omega;
    M[50,58] = omega;
    M[51,6] = omega;
    M[51,19] = omega;
    M[51,35] = omega * kappa;
    M[51,47] = omega * kappa;
    M[51,50] = kappa;
    M[51,52] = 1;
    M[51,53] = 1;
    M[51,55] = omega;
    M[51,59] = omega;
    M[52,7] = omega;
    M[52,20] = omega;
    M[52,36] = omega * kappa;
    M[52,48] = omega * kappa;
    M[52,50] = 1;
    M[52,51] = 1;
    M[52,53] = kappa;
    M[52,56] = omega;
    M[52,60] = omega;
    M[53,8] = omega;
    M[53,21] = omega;
    M[53,37] = omega * kappa;
    M[53,49] = omega * kappa;
    M[53,50] = 1;
    M[53,51] = 1;
    M[53,52] = kappa;
    M[53,57] = omega;
    M[53,61] = omega;
    M[54,9] = omega;
    M[54,22] = omega;
    M[54,38] = omega * kappa;
    M[54,46] = omega;
    M[54,50] = omega;
    M[54,55] = kappa;
    M[54,56] = omega;
    M[54,57] = omega;
    M[54,58] = omega * kappa;
    M[55,10] = omega;
    M[55,23] = omega;
    M[55,39] = omega * kappa;
    M[55,47] = omega;
    M[55,51] = omega;
    M[55,54] = kappa;
    M[55,56] = omega;
    M[55,57] = omega;
    M[55,59] = omega * kappa;
    M[56,24] = omega;
    M[56,40] = omega * kappa;
    M[56,48] = omega;
    M[56,52] = omega;
    M[56,54] = omega;
    M[56,55] = omega;
    M[56,57] = kappa;
    M[56,60] = omega * kappa;
    M[57,25] = omega;
    M[57,41] = omega * kappa;
    M[57,49] = omega;
    M[57,53] = omega;
    M[57,54] = omega;
    M[57,55] = omega;
    M[57,56] = kappa;
    M[57,61] = omega * kappa;
    M[58,11] = omega;
    M[58,26] = omega;
    M[58,42] = omega * kappa;
    M[58,46] = omega;
    M[58,50] = omega;
    M[58,54] = omega * kappa;
    M[58,59] = kappa;
    M[58,60] = 1;
    M[58,61] = 1;
    M[59,12] = omega;
    M[59,27] = omega;
    M[59,43] = omega * kappa;
    M[59,47] = omega;
    M[59,51] = omega;
    M[59,55] = omega * kappa;
    M[59,58] = kappa;
    M[59,60] = 1;
    M[59,61] = 1;
    M[60,28] = omega;
    M[60,44] = omega * kappa;
    M[60,48] = omega;
    M[60,52] = omega;
    M[60,56] = omega * kappa;
    M[60,58] = 1;
    M[60,59] = 1;
    M[60,61] = kappa;
    M[61,13] = omega;
    M[61,29] = omega;
    M[61,45] = omega * kappa;
    M[61,49] = omega;
    M[61,53] = omega;
    M[61,57] = omega * kappa;
    M[61,58] = 1;
    M[61,59] = 1;
    M[61,60] = kappa;
    
    
    // Fill in lower diagonal
    // M += M'; 
    
    // Scale by mu
    // M *= mu;
    
    // Multiply by equilibrium frequencies
    // *Apply the equilibrium frequencies*/
    for(i in 1:61) {
      sqp = sqrt(pi_eq[i]);
      for(j in 1:61) {
        M[i][j] *= sqp;
        M[j][i] *= sqp;
      }
    }
    
    for(i in 1:61){
      for(j in 1:61){ 
        M[i, i] -= M[i, j]*sqrt(pi_eq[j] * pi_eq[i]);
      }
    }
    
    // Compute the diagonal
    // M = add_diag(M, -(M * rep_vector(1, 61)));
    
    return(M);
  }
  
  real partial_sum_constant(array[] int X_slice,
  int start, int end,
  matrix m_AB,
  array[,] int X) {
    
    vector[61] logvec;
    vector[end - start + 1] lv;
    int z = 1; 
    real sqp, mABs, n, ll;
    row_vector[61] mAB;
    row_vector[61] obs;
    
    for(loc in start:end){
      // Likelihood calculation
      profile("likelihood"){
        obs = to_row_vector(X[loc,]);
        n = sum(obs);
        for(i in 1:61){
          mAB = m_AB[i, ]; 
          mABs = sum(m_AB[i, ]);
          // Dirichlet-multinomial lpmf for this ancestral codon
          ll = - (lgamma(mABs + n) - lgamma(mABs));
          for(k in 1:61){
            ll += (lgamma(mAB[k] + obs[k]) - lgamma(mAB[k]));
          }
          logvec[i] = ll;
        }
      }
      // Sum likelihood over all ancestors
      lv[z] = sum(logvec);
      z += 1;
    }
    
    return sum(lv);
  }
  
  real partial_sum(array[] int X_slice,
  int start, int end,
  row_vector pi_eq,
  vector mu,
  vector omega,
  vector kappa,
  array[,] int X) {
    
    real m_AA;
    matrix[61, 61] m_AB;
    matrix[61, 61] mutmat;
    matrix[61, 61] V;
    matrix[61, 61] Va;
    vector[61] D;
    matrix[61, 61] V_inv;
    vector[61] logvec;
    vector[end - start + 1] lv;
    int z = 1; 
    real sqp;
    real mABs, sobs;
    row_vector[61] obs;
    
    for(loc in start:end){
      // Calculate substitution rate matrix
      profile("mutmat") {
        mutmat = PDRM(mu[loc], kappa[loc], omega[loc], pi_eq);
      }
      
      // Eigenvectors/values of substitution rate matrix
      profile("eigen") {
        V = eigenvectors_sym(mutmat);
        D = 1 / (1 - eigenvalues_sym(mutmat));
        V_inv = diag_post_multiply(V, D);
      }
      
      profile("m_AB"){
        // Create m_AB for each ancestral codon
        for(i in 1:61) {
          Va = rep_matrix(row(V, i), 61);
          m_AB[i, ] = to_row_vector(rows_dot_product(Va, V_inv));
        }
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
        }
      }
      
      // Likelihood calculation
      profile("likelihood"){
        obs = to_row_vector(X[loc,]);
        sobs = sum(obs);
        for(i in 1:61){
          mABs = sum(m_AB[i, ]);
          logvec[i] = (lgamma(mABs) + sum(lgamma(m_AB[i, ] + obs)) - lgamma(mABs + sobs) - sum(lgamma(m_AB[i, ])));
        }
      }
      // Sum likelihood over all ancestors
      lv[z] = sum(logvec);
      z += 1;
    }
    
    return sum(lv);
  }
  
  // substitution rate matrix
  matrix PDRM(real mu, real kappa, real omega, row_vector pi_eq) {
    matrix[61,61] M = rep_matrix(0.,61,61);
    real sqp;
    
    M[1,2] = kappa;
    M[1,3] = omega;
    M[1,4] = omega;
    M[1,5] = omega * kappa;
    M[1,9] = omega;
    M[1,11] = omega;
    M[1,14] = omega * kappa;
    M[1,30] = omega;
    M[1,46] = omega;
    M[2,1] = kappa;
    M[2,3] = omega;
    M[2,4] = omega;
    M[2,6] = omega * kappa;
    M[2,10] = omega;
    M[2,12] = omega;
    M[2,15] = omega * kappa;
    M[2,31] = omega;
    M[2,47] = omega;
    M[3,1] = omega;
    M[3,2] = omega;
    M[3,4] = kappa;
    M[3,7] = omega * kappa;
    M[3,16] = kappa;
    M[3,32] = omega;
    M[3,48] = omega;
    M[4,1] = omega;
    M[4,2] = omega;
    M[4,3] = kappa;
    M[4,8] = omega * kappa;
    M[4,13] = omega;
    M[4,17] = kappa;
    M[4,33] = omega;
    M[4,49] = omega;
    M[5,1] = omega * kappa;
    M[5,6] = kappa;
    M[5,7] = 1;
    M[5,8] = 1;
    M[5,9] = omega;
    M[5,11] = omega;
    M[5,18] = omega * kappa;
    M[5,34] = omega;
    M[5,50] = omega;
    M[6,2] = omega * kappa;
    M[6,5] = kappa;
    M[6,7] = 1;
    M[6,8] = 1;
    M[6,10] = omega;
    M[6,12] = omega;
    M[6,19] = omega * kappa;
    M[6,35] = omega;
    M[6,51] = omega;
    M[7,3] = omega * kappa;
    M[7,5] = 1;
    M[7,6] = 1;
    M[7,8] = kappa;
    M[7,20] = omega * kappa;
    M[7,36] = omega;
    M[7,52] = omega;
    M[8,4] = omega * kappa;
    M[8,5] = 1;
    M[8,6] = 1;
    M[8,7] = kappa;
    M[8,13] = omega;
    M[8,21] = omega * kappa;
    M[8,37] = omega;
    M[8,53] = omega;
    M[9,1] = omega;
    M[9,5] = omega;
    M[9,10] = kappa;
    M[9,11] = omega * kappa;
    M[9,22] = omega * kappa;
    M[9,38] = omega;
    M[9,54] = omega;
    M[10,2] = omega;
    M[10,6] = omega;
    M[10,9] = kappa;
    M[10,12] = omega * kappa;
    M[10,23] = omega * kappa;
    M[10,39] = omega;
    M[10,55] = omega;
    M[11,1] = omega;
    M[11,5] = omega;
    M[11,9] = omega * kappa;
    M[11,12] = kappa;
    M[11,13] = omega;
    M[11,26] = omega * kappa;
    M[11,42] = omega;
    M[11,58] = omega;
    M[12,2] = omega;
    M[12,6] = omega;
    M[12,10] = omega * kappa;
    M[12,11] = kappa;
    M[12,13] = omega;
    M[12,27] = omega * kappa;
    M[12,43] = omega;
    M[12,59] = omega;
    M[13,4] = omega;
    M[13,8] = omega;
    M[13,11] = omega;
    M[13,12] = omega;
    M[13,29] = omega * kappa;
    M[13,45] = omega;
    M[13,61] = omega;
    M[14,1] = omega * kappa;
    M[14,15] = kappa;
    M[14,16] = 1;
    M[14,17] = 1;
    M[14,18] = omega * kappa;
    M[14,22] = omega;
    M[14,26] = omega;
    M[14,30] = omega;
    M[14,46] = omega;
    M[15,2] = omega * kappa;
    M[15,14] = kappa;
    M[15,16] = 1;
    M[15,17] = 1;
    M[15,19] = omega * kappa;
    M[15,23] = omega;
    M[15,27] = omega;
    M[15,31] = omega;
    M[15,47] = omega;
    M[16,3] = kappa;
    M[16,14] = 1;
    M[16,15] = 1;
    M[16,17] = kappa;
    M[16,20] = omega * kappa;
    M[16,24] = omega;
    M[16,28] = omega;
    M[16,32] = omega;
    M[16,48] = omega;
    M[17,4] = kappa;
    M[17,14] = 1;
    M[17,15] = 1;
    M[17,16] = kappa;
    M[17,21] = omega * kappa;
    M[17,25] = omega;
    M[17,29] = omega;
    M[17,33] = omega;
    M[17,49] = omega;
    M[18,5] = omega * kappa;
    M[18,14] = omega * kappa;
    M[18,19] = kappa;
    M[18,20] = 1;
    M[18,21] = 1;
    M[18,22] = omega;
    M[18,26] = omega;
    M[18,34] = omega;
    M[18,50] = omega;
    M[19,6] = omega * kappa;
    M[19,15] = omega * kappa;
    M[19,18] = kappa;
    M[19,20] = 1;
    M[19,21] = 1;
    M[19,23] = omega;
    M[19,27] = omega;
    M[19,35] = omega;
    M[19,51] = omega;
    M[20,7] = omega * kappa;
    M[20,16] = omega * kappa;
    M[20,18] = 1;
    M[20,19] = 1;
    M[20,21] = kappa;
    M[20,24] = omega;
    M[20,28] = omega;
    M[20,36] = omega;
    M[20,52] = omega;
    M[21,8] = omega * kappa;
    M[21,17] = omega * kappa;
    M[21,18] = 1;
    M[21,19] = 1;
    M[21,20] = kappa;
    M[21,25] = omega;
    M[21,29] = omega;
    M[21,37] = omega;
    M[21,53] = omega;
    M[22,9] = omega * kappa;
    M[22,14] = omega;
    M[22,18] = omega;
    M[22,23] = kappa;
    M[22,24] = omega;
    M[22,25] = omega;
    M[22,26] = omega * kappa;
    M[22,38] = omega;
    M[22,54] = omega;
    M[23,10] = omega * kappa;
    M[23,15] = omega;
    M[23,19] = omega;
    M[23,22] = kappa;
    M[23,24] = omega;
    M[23,25] = omega;
    M[23,27] = omega * kappa;
    M[23,39] = omega;
    M[23,55] = omega;
    M[24,16] = omega;
    M[24,20] = omega;
    M[24,22] = omega;
    M[24,23] = omega;
    M[24,25] = kappa;
    M[24,28] = omega * kappa;
    M[24,40] = omega;
    M[24,56] = omega;
    M[25,17] = omega;
    M[25,21] = omega;
    M[25,22] = omega;
    M[25,23] = omega;
    M[25,24] = kappa;
    M[25,29] = omega * kappa;
    M[25,41] = omega;
    M[25,57] = omega;
    M[26,11] = omega * kappa;
    M[26,14] = omega;
    M[26,18] = omega;
    M[26,22] = omega * kappa;
    M[26,27] = kappa;
    M[26,28] = 1;
    M[26,29] = 1;
    M[26,42] = omega;
    M[26,58] = omega;
    M[27,12] = omega * kappa;
    M[27,15] = omega;
    M[27,19] = omega;
    M[27,23] = omega * kappa;
    M[27,26] = kappa;
    M[27,28] = 1;
    M[27,29] = 1;
    M[27,43] = omega;
    M[27,59] = omega;
    M[28,16] = omega;
    M[28,20] = omega;
    M[28,24] = omega * kappa;
    M[28,26] = 1;
    M[28,27] = 1;
    M[28,29] = kappa;
    M[28,44] = 1;
    M[28,60] = omega;
    M[29,13] = omega * kappa;
    M[29,17] = omega;
    M[29,21] = omega;
    M[29,25] = omega * kappa;
    M[29,26] = 1;
    M[29,27] = 1;
    M[29,28] = kappa;
    M[29,45] = 1;
    M[29,61] = omega;
    M[30,1] = omega;
    M[30,14] = omega;
    M[30,31] = kappa;
    M[30,32] = 1;
    M[30,33] = omega;
    M[30,34] = omega * kappa;
    M[30,38] = omega;
    M[30,42] = omega;
    M[30,46] = omega * kappa;
    M[31,2] = omega;
    M[31,15] = omega;
    M[31,30] = kappa;
    M[31,32] = 1;
    M[31,33] = omega;
    M[31,35] = omega * kappa;
    M[31,39] = omega;
    M[31,43] = omega;
    M[31,47] = omega * kappa;
    M[32,3] = omega;
    M[32,16] = omega;
    M[32,30] = 1;
    M[32,31] = 1;
    M[32,33] = omega * kappa;
    M[32,36] = omega * kappa;
    M[32,40] = omega;
    M[32,44] = omega;
    M[32,48] = omega * kappa;
    M[33,4] = omega;
    M[33,17] = omega;
    M[33,30] = omega;
    M[33,31] = omega;
    M[33,32] = omega * kappa;
    M[33,37] = omega * kappa;
    M[33,41] = omega;
    M[33,45] = omega;
    M[33,49] = omega * kappa;
    M[34,5] = omega;
    M[34,18] = omega;
    M[34,30] = omega * kappa;
    M[34,35] = kappa;
    M[34,36] = 1;
    M[34,37] = 1;
    M[34,38] = omega;
    M[34,42] = omega;
    M[34,50] = omega * kappa;
    M[35,6] = omega;
    M[35,19] = omega;
    M[35,31] = omega * kappa;
    M[35,34] = kappa;
    M[35,36] = 1;
    M[35,37] = 1;
    M[35,39] = omega;
    M[35,43] = omega;
    M[35,51] = omega * kappa;
    M[36,7] = omega;
    M[36,20] = omega;
    M[36,32] = omega * kappa;
    M[36,34] = 1;
    M[36,35] = 1;
    M[36,37] = kappa;
    M[36,40] = omega;
    M[36,44] = omega;
    M[36,52] = omega * kappa;
    M[37,8] = omega;
    M[37,21] = omega;
    M[37,33] = omega * kappa;
    M[37,34] = 1;
    M[37,35] = 1;
    M[37,36] = kappa;
    M[37,41] = omega;
    M[37,45] = omega;
    M[37,53] = omega * kappa;
    M[38,9] = omega;
    M[38,22] = omega;
    M[38,30] = omega;
    M[38,34] = omega;
    M[38,39] = kappa;
    M[38,40] = omega;
    M[38,41] = omega;
    M[38,42] = omega * kappa;
    M[38,54] = omega * kappa;
    M[39,10] = omega;
    M[39,23] = omega;
    M[39,31] = omega;
    M[39,35] = omega;
    M[39,38] = kappa;
    M[39,40] = omega;
    M[39,41] = omega;
    M[39,43] = omega * kappa;
    M[39,55] = omega * kappa;
    M[40,24] = omega;
    M[40,32] = omega;
    M[40,36] = omega;
    M[40,38] = omega;
    M[40,39] = omega;
    M[40,41] = kappa;
    M[40,44] = omega * kappa;
    M[40,56] = omega * kappa;
    M[41,25] = omega;
    M[41,33] = omega;
    M[41,37] = omega;
    M[41,38] = omega;
    M[41,39] = omega;
    M[41,40] = kappa;
    M[41,45] = omega * kappa;
    M[41,57] = omega * kappa;
    M[42,11] = omega;
    M[42,26] = omega;
    M[42,30] = omega;
    M[42,34] = omega;
    M[42,38] = omega * kappa;
    M[42,43] = kappa;
    M[42,44] = omega;
    M[42,45] = omega;
    M[42,58] = omega * kappa;
    M[43,12] = omega;
    M[43,27] = omega;
    M[43,31] = omega;
    M[43,35] = omega;
    M[43,39] = omega * kappa;
    M[43,42] = kappa;
    M[43,44] = omega;
    M[43,45] = omega;
    M[43,59] = omega * kappa;
    M[44,28] = 1;
    M[44,32] = omega;
    M[44,36] = omega;
    M[44,40] = omega * kappa;
    M[44,42] = omega;
    M[44,43] = omega;
    M[44,45] = kappa;
    M[44,60] = omega * kappa;
    M[45,13] = omega;
    M[45,29] = 1;
    M[45,33] = omega;
    M[45,37] = omega;
    M[45,41] = omega * kappa;
    M[45,42] = omega;
    M[45,43] = omega;
    M[45,44] = kappa;
    M[45,61] = omega * kappa;
    M[46,1] = omega;
    M[46,14] = omega;
    M[46,30] = omega * kappa;
    M[46,47] = kappa;
    M[46,48] = 1;
    M[46,49] = 1;
    M[46,50] = omega * kappa;
    M[46,54] = omega;
    M[46,58] = omega;
    M[47,2] = omega;
    M[47,15] = omega;
    M[47,31] = omega * kappa;
    M[47,46] = kappa;
    M[47,48] = 1;
    M[47,49] = 1;
    M[47,51] = omega * kappa;
    M[47,55] = omega;
    M[47,59] = omega;
    M[48,3] = omega;
    M[48,16] = omega;
    M[48,32] = omega * kappa;
    M[48,46] = 1;
    M[48,47] = 1;
    M[48,49] = kappa;
    M[48,52] = omega * kappa;
    M[48,56] = omega;
    M[48,60] = omega;
    M[49,4] = omega;
    M[49,17] = omega;
    M[49,33] = omega * kappa;
    M[49,46] = 1;
    M[49,47] = 1;
    M[49,48] = kappa;
    M[49,53] = omega * kappa;
    M[49,57] = omega;
    M[49,61] = omega;
    M[50,5] = omega;
    M[50,18] = omega;
    M[50,34] = omega * kappa;
    M[50,46] = omega * kappa;
    M[50,51] = kappa;
    M[50,52] = 1;
    M[50,53] = 1;
    M[50,54] = omega;
    M[50,58] = omega;
    M[51,6] = omega;
    M[51,19] = omega;
    M[51,35] = omega * kappa;
    M[51,47] = omega * kappa;
    M[51,50] = kappa;
    M[51,52] = 1;
    M[51,53] = 1;
    M[51,55] = omega;
    M[51,59] = omega;
    M[52,7] = omega;
    M[52,20] = omega;
    M[52,36] = omega * kappa;
    M[52,48] = omega * kappa;
    M[52,50] = 1;
    M[52,51] = 1;
    M[52,53] = kappa;
    M[52,56] = omega;
    M[52,60] = omega;
    M[53,8] = omega;
    M[53,21] = omega;
    M[53,37] = omega * kappa;
    M[53,49] = omega * kappa;
    M[53,50] = 1;
    M[53,51] = 1;
    M[53,52] = kappa;
    M[53,57] = omega;
    M[53,61] = omega;
    M[54,9] = omega;
    M[54,22] = omega;
    M[54,38] = omega * kappa;
    M[54,46] = omega;
    M[54,50] = omega;
    M[54,55] = kappa;
    M[54,56] = omega;
    M[54,57] = omega;
    M[54,58] = omega * kappa;
    M[55,10] = omega;
    M[55,23] = omega;
    M[55,39] = omega * kappa;
    M[55,47] = omega;
    M[55,51] = omega;
    M[55,54] = kappa;
    M[55,56] = omega;
    M[55,57] = omega;
    M[55,59] = omega * kappa;
    M[56,24] = omega;
    M[56,40] = omega * kappa;
    M[56,48] = omega;
    M[56,52] = omega;
    M[56,54] = omega;
    M[56,55] = omega;
    M[56,57] = kappa;
    M[56,60] = omega * kappa;
    M[57,25] = omega;
    M[57,41] = omega * kappa;
    M[57,49] = omega;
    M[57,53] = omega;
    M[57,54] = omega;
    M[57,55] = omega;
    M[57,56] = kappa;
    M[57,61] = omega * kappa;
    M[58,11] = omega;
    M[58,26] = omega;
    M[58,42] = omega * kappa;
    M[58,46] = omega;
    M[58,50] = omega;
    M[58,54] = omega * kappa;
    M[58,59] = kappa;
    M[58,60] = 1;
    M[58,61] = 1;
    M[59,12] = omega;
    M[59,27] = omega;
    M[59,43] = omega * kappa;
    M[59,47] = omega;
    M[59,51] = omega;
    M[59,55] = omega * kappa;
    M[59,58] = kappa;
    M[59,60] = 1;
    M[59,61] = 1;
    M[60,28] = omega;
    M[60,44] = omega * kappa;
    M[60,48] = omega;
    M[60,52] = omega;
    M[60,56] = omega * kappa;
    M[60,58] = 1;
    M[60,59] = 1;
    M[60,61] = kappa;
    M[61,13] = omega;
    M[61,29] = omega;
    M[61,45] = omega * kappa;
    M[61,49] = omega;
    M[61,53] = omega;
    M[61,57] = omega * kappa;
    M[61,58] = 1;
    M[61,59] = 1;
    M[61,60] = kappa;
    
    
    // Fill in lower diagonal
    // M += M'; 
    
    // Scale by mu
    M *= mu;
    
    // Multiply by equilibrium frequencies
    for(i in 1:61){
      M[i, ] .*= pi_eq;
    }
    
    // Compute the diagonal
    M = add_diag(M, -(M * rep_vector(1, 61)));
    
    return(M);
  }
  
    // real dirichlet_multinomial_lpmf(array[] int y, vector alpha) {
  //   real alpha_plus = sum(alpha);
  //   return lgamma(alpha_plus) + sum(lgamma(alpha + to_vector(y)))
  //   - lgamma(alpha_plus+sum(y)) - sum(lgamma(alpha));
  // }
  
  
  // real partial_sum_wild(array[] int X_slice,
  // int start, int end,
  // row_vector pi_eq,
  // vector mu,
  // vector omega,
  // vector kappa,
  // array[,] int X,
  // array[,] int A,
  // array[] int csac,
  // array[] int nac) {
    //   
    //   vector[61] m_AA;
    //   vector[61] m_AB;
    //   vector[61] alpha_A;
    //   matrix[61,61] mutmat;
    //   matrix[61,61] V;
    //   vector[61] D;
    //   matrix[61, 61] V_inv;
    //   array[61] matrix[61,61] Va;
    //   vector[61] logvec = rep_vector(0, 61);
    //   vector[end - start + 1] lv = rep_vector(0, end - start + 1);
    //   int z = 1;
    //   int cur_anc; 
    //   
    //   for(i in start:end){
      //     // Calculate substitution rate matrix
      //     profile("mutmat") {
        //       mutmat = PDRM(mu[i], kappa[i], omega[i], pi_eq);
        //     }
        //     // Eigenvectors/values of substitution rate matrix
        //     profile("eigen") {
          //       V = eigenvectors_sym(mutmat);
          //       D = 1 / (1 - eigenvalues_sym(mutmat));
          //       // V^-1 = transpose(V) when V symmetric
          //       V_inv = diag_post_multiply(V, D)';
          //       
          //       for(ii in 1:61){
            //         Va[ii] = rep_matrix(row(V, ii), 61);
            //       }
            //     }
            //     
            //     m_AA = rows_dot_product(V, V_inv);
            //     // Loop over ancestral codons
            //     for(k in 1:nac[i]){
              //       profile("m_AB") {
                //         // Calculates a vector of m_AB where A is the ancestral codon (k)
                //         // and the index of the vector is the codon B
                //         // Normalise vector by m_AA
                //         cur_anc = A[(csac[i] + (k - 1)), 2];
                //         m_AB = rows_dot_product(Va[cur_anc], V_inv) / m_AA[cur_anc];
                //         
                //       }
                //       profile("likelihood") {
                  //         logvec[k] = dirichlet_multinomial_lpmf(X[i, ] | m_AB);
                  //       }
                  //     }
                  //     
                  //     // Vector allocation then sum() is apparently more efficient than incrementing within a loop
                  //     lv[z] = sum(logvec);
                  //     logvec = rep_vector(0, 61);
                  //     z += 1;
                  //   }
                  // 
                  //   return sum(lv);
                  // }
                  
}
