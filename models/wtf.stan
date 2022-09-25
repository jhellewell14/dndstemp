functions {
  #include functions/GTR_mat.stan
  #include functions/calc_mAB.stan
}

data {
  row_vector[61] pi_eq;
  int <lower = 1> l; // gene length
  array[l, 61] int X; // codon frequencies at each locus in the gene
  array [l] int N;
  int <lower = 1> K; // number of mixture distributions
  vector[K] om_priors;
  int <lower = 1> M; // number of shards (cores)
  array[M] int n_per_shard;
  int max_per_shard;
  array[M] int shard_starts;
  array[M] int shard_ends;
  vector[61 * 61] map_pars;
} 


transformed data{
  row_vector[K] omega = to_row_vector(om_priors);
  vector[61] lp = to_vector(log(pi_eq));
  array[l] matrix[61, 61] obs_mat;
  matrix[61, 61] pimat = diag_matrix(sqrt(to_vector(pi_eq)));
  matrix[61, 61] pimatinv = diag_matrix(inv(sqrt(to_vector(pi_eq))));
  matrix[61, 61] pimult;
  array[l] real phi;
  array[61] real pi_array = to_array_1d(lp);
  vector[61] ones = rep_vector(1, 61);
  array[M, ((max_per_shard * 61) + 4)] int xi;
  array [M, 61 + l] real xr;
  array[M] vector[0] thetay;
  
  for(j in 1:61){
    for(i in 1:61){
      pimult[i, j] = sqrt(pi_eq[j] / pi_eq[i]);
    }
  }
  
  for(i in 1:l){
    phi[i] = lgamma(N[i] + 1) - sum(lgamma(to_vector(X[i, ]) + 1));
    obs_mat[i] = rep_matrix(to_row_vector(X[i, ]), 61);
  }
  
  
  int st, en;
  for(i in 1:M){
    // Fill xi
    for(j in 1:max_per_shard){
      st = 1 + ((j - 1) * 61);
      en = st + 60;
      if(j <= n_per_shard[i]){
        xi[i, st:en] = X[shard_starts[i] + (j - 1),];
      }else{
        xi[i ,st:en] = rep_array(-1, 61);
      }
      
    }
    xi[i, (max_per_shard * 61) + 1] = n_per_shard[i];
    xi[i, (max_per_shard * 61) + 2] = K;
    xi[i, (max_per_shard * 61) + 3] = shard_starts[i];
    xi[i, (max_per_shard * 61) + 4] = shard_ends[i];
    
    // Fill xr
    xr[i, ] = append_array(pi_array, phi);
  }
}

generated quantities {
  vector[l] log_lik = map_rect(f2, map_pars, thetay, xr, xi);
  // print(log_lik);
  // vector[n_per_shard[1]] l1 = f2(map_pars, thetay[1], xr[1], xi[1]);
  // print(l1);
  // print("hello");
  // vector[n_per_shard[2]] l2;
  // vector[n_per_shard[3]] l3;
  // vector[n_per_shard[4]] l4;
  // vector[n_per_shard[5]] l5;
  // vector[n_per_shard[6]] l6;
  // vector[n_per_shard[7]] l7;
  // vector[n_per_shard[8]] l8;
  // 
  // profile("ll"){
  //   l1 = f2(map_pars, thetay[1], xr[1], xi[1]);
  //   l2 = f2(map_pars, thetay[1], xr[2], xi[2]);
  //   l3 = f2(map_pars, thetay[1], xr[3], xi[3]);
  //   l4 = f2(map_pars, thetay[1], xr[4], xi[4]);
  //   l5 = f2(map_pars, thetay[1], xr[5], xi[5]);
  //   l6 = f2(map_pars, thetay[1], xr[6], xi[6]);
  //   l7 = f2(map_pars, thetay[1], xr[7], xi[7]);
  //   l8 = f2(map_pars, thetay[1], xr[8], xi[8]);
  // }
}
