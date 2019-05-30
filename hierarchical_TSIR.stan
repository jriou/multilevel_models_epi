data {
  int inference;
  int<lower=1> W; // number of records
  int<lower=1> K; // number of islands
  int<lower=1> J; // number of viruses
  int<lower=1> island[W]; // identification of island 1=TAHITI 2=TUAMOTU 3=MOOREA 4=SLV 5=MARQUISES 6=AUSTRALES 7=SAINT MARTIN 8=MARTINIQUE 9=GUADELOUPE
  int<lower=0> virus[W]; // identification of virus  0=CHIKV 1=ZIKV
  int<lower=0> O_t[W]; // number of reported cases at time t
  real<lower=0> Ostar_t[W]; // exposure at time t
  int<lower=0> sumO_t[W]; // cumulative number of reported cases at time t
  int<lower=0> sumO_tmaxz[K]; // total number of reported ZIKV cases
  int<lower=0> sumO_tmaxc[K]; // total number of reported CHIKV cases
  int<lower=0> pop[K]; // island population
  int<lower=0> C; // number of weather covariables
  matrix[W,C] weather; // matrix of weather covariables 
}

parameters {
// hyperparameters
  real mu_b0c; // distribution of b0 for CHIKV
  real<lower=0> sigma_b0c;
  real mu_rhoc; // logit of the distribution of rho for CHIKV
  real<lower=0> sigma_rhoc;

// random parameters
  real b0ic[K]; // island-specific base transmission for CHIKV
  real logitrhoc[K]; // island-specific logit reporting rate for CHIKV

// fixed parameters
  real bD; // difference of transmission for ZIKV
  real omega; // difference of reporting for ZIKV
  vector[C] bW; // effects of the weather covariables

// dispersion par
  real<lower=0,upper=100> phi;
}

transformed parameters {
  real<lower=0,upper=1> rhoc[K]; // reporting rate
  real<lower=0,upper=1> rhoz[K]; // reporting rate by island
  vector[W] theta;
  real<lower=0> beta[W];
  real<lower=0> lp[W];
  real sampledisp[W];

  for (i in 1:K) {
    rhoc[i] = exp(logitrhoc[i]) / (1+exp(logitrhoc[i]));
    rhoz[i] = rhoc[i]*exp(omega)/(1-rhoc[i]+rhoc[i]*exp(omega));
  }
  theta = weather * bW;
  for( i in 1:W) {
    beta[i] = exp( b0ic[island[i]] + bD*virus[i] + theta[i] );
    lp[i] = beta[i] * Ostar_t[i] * ( 1 - sumO_t[i] / ( if_else(virus[i]==1,rhoz[island[i]],rhoc[island[i]]) * pop[island[i]]) );
    sampledisp[i] = lp[i]/phi;
  }
}


model {
// hyperparameter priors 
  mu_b0c ~ student_t(5,0,2.5);
  sigma_b0c ~ cauchy(0,2.5);
  mu_rhoc ~ normal(0,1);
  sigma_rhoc ~ inv_gamma(10,10);

// random parameter priors
  for (i in 1:K) {
    b0ic[i] ~ normal(mu_b0c,sigma_b0c);
    logitrhoc[i] ~ normal(mu_rhoc,sigma_rhoc);
  }

// fixed parameters priors
  bD ~ student_t(5,0,2.5); 
  omega ~ normal(0,1);
  for (i in 1:C) {
    bW[i] ~ student_t(5,0,2.5);
  }

// likelihood
  if(inference==1) target += neg_binomial_2_lpmf(O_t|lp,sampledisp);
}


generated quantities {
  real log_lik[W]; // save the likelihood for LOOIC computation
  for (i in 1:W) {
    log_lik[i] = neg_binomial_2_lpmf(O_t[i]|lp[i],sampledisp[i]);
  }
}
