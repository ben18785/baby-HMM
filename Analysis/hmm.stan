data {
  int<lower=0> N;
  int<lower=0,upper=120> simple[N];
  int<lower=0,upper=120> complex[N];
  int<lower=1> K;
}

transformed data{
 matrix[3,3] beta;
 vector[3] zeta_1;
 vector[3] zeta_2;
 beta[1,1] = 10;
 beta[1,2] = 1;
 beta[1,3] = 1;
 beta[2,1] = 1;
 beta[2,2] = 10;
 beta[2,3] = 1;
 beta[3,1] = 1;
 beta[3,2] = 1;
 beta[3,3] = 10;
 zeta_1[1] = 10;
 zeta_1[2] = 10;
 zeta_1[3] = 110;
 zeta_2[1] = 10;
 zeta_2[2] = 110;
 zeta_2[3] = 10;
}

parameters {
  vector[K] mu1;
  real<lower=0> sigma1[K];
  vector[K] muTemp;
  real<lower=0> sigma2[K];
  simplex[K] theta[K];
}

transformed parameters {
  vector[K] log_theta_tr[K];
  real mu2[K];
  
  for(i in 1:K)
    mu2[i] = 120 - muTemp[i];
  
   for (k_from in 1:K)
    for (k in 1:K)
      log_theta_tr[k, k_from] = log(theta[k_from, k]);
    
}

model {
  vector[K] lp;
  vector[K] lp_p1;
  
  lp = rep_vector(-log(K), K);
  
  // Forwards algorithm
  for (n in 1:N) {
    real logSimple[K];
      real logComplex[K];
      real logTheta[K];
    for (k in 1:K){
      logTheta[k] = log_sum_exp(log_theta_tr[k] + lp);
      if((simple[n] > 0) && (simple[n] < 120)){
        logSimple[k] = normal_lpdf(simple[n]|mu1[k],sigma1[k]); 
      } else if(simple[n] < 0.0001){
        logSimple[k] = normal_lcdf(0|mu1[k],sigma1[k]);
      } else if(simple[n] > 119.9999){
        logSimple[k] = normal_lccdf(120|mu1[k],sigma1[k]);
      }
      if((complex[n] > 0) && (complex[n] < 120)){
        logComplex[k] = normal_lpdf(complex[n]|120-muTemp[k],sigma2[k]); 
      } else if(complex[n] < 0.0001){
        logComplex[k] = normal_lcdf(0|120-muTemp[k],sigma2[k]);
      } else if(complex[n] > 119.9999){
        logComplex[k] = normal_lccdf(120|120-muTemp[k],sigma2[k]);
      }
      if(n == 1){
        lp_p1[k] = logSimple[k] + logComplex[k];
      }else{
        lp_p1[k] = logTheta[k] + logSimple[k] + logComplex[k];
      }
    }
    lp = lp_p1;
  }
  target += log_sum_exp(lp);
  for(i in 1:K){
    mu1[i] ~ normal(zeta_1[i],20);
    mu2[i] ~ normal(zeta_2[i],20);
  }
    
  sigma1 ~ normal(20,5);
  sigma2 ~ normal(20,5);
  
  // for(i in 1:K)
  //   theta[i] ~ dirichlet(to_vector(beta[i]));
}

generated quantities{

  //// Viterbi to get states
  int<lower=1,upper=K> state[N];
  real log_p_y_star;
  {
    int back_ptr[N, K];
    real best_logp[N, K];
    real best_total_logp;
    real logSimple[K];
      real logComplex[K];
      real logTheta[K];
    for (k in 1:K){
      if((simple[1] > 0) && (simple[1] < 120)){
        logSimple[k] = normal_lpdf(simple[1]|mu1[k],sigma1[k]); 
      } else if(simple[1] < 0.0001){
        logSimple[k] = normal_lcdf(0|mu1[k],sigma1[k]);
      } else if(simple[1] > 119.9999){
        logSimple[k] = normal_lccdf(120|mu1[k],sigma1[k]);
      }
      if((complex[1] > 0) && (complex[1] < 120)){
        logComplex[k] = normal_lpdf(complex[1]|120-muTemp[k],sigma2[k]); 
      } else if(complex[1] < 0.0001){
        logComplex[k] = normal_lcdf(0|120-muTemp[k],sigma2[k]);
      } else if(complex[1] > 119.9999){
        logComplex[k] = normal_lccdf(120|120-muTemp[k],sigma2[k]);
      }
      best_logp[1, k] = logSimple[k] + logComplex[k];
    }
    for (t in 2:N) {
      for (k in 1:K) {
      best_logp[t, k] = negative_infinity();
        for (j in 1:K) {
          real logp;
          
      if((simple[t] > 0) && (simple[t] < 120)){
        logSimple[k] = normal_lpdf(simple[1]|mu1[k],sigma1[k]); 
      } else if(simple[t] < 0.0001){
        logSimple[k] = normal_lcdf(0|mu1[k],sigma1[k]);
      } else if(simple[t] > 119.9999){
        logSimple[k] = normal_lccdf(120|mu1[k],sigma1[k]);
      }
      if((complex[t] > 0) && (complex[t] < 120)){
        logComplex[k] = normal_lpdf(complex[t]|120-muTemp[k],sigma2[k]); 
      } else if(complex[t] < 0.0001){
        logComplex[k] = normal_lcdf(0|120-muTemp[k],sigma2[k]);
      } else if(complex[t] > 119.9999){
        logComplex[k] = normal_lccdf(120|120-muTemp[k],sigma2[k]);
      }
      logp = best_logp[t-1, j] + logSimple[k] + logComplex[k] + log(theta[j,k]);
          if (logp > best_logp[t, k]) {
            back_ptr[t, k] = j;
            best_logp[t, k] = logp;
          }
        }
      }
    }
    log_p_y_star = max(best_logp[N]);
    for (k in 1:K)
      if (best_logp[N, k] == log_p_y_star)
      state[N] = k;
      for (t in 1:(N - 1))
      state[N - t] = back_ptr[N - t + 1,state[N - t + 1]];
  }

}