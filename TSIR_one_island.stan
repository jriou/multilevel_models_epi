data {
  int inference;
  int<lower=1> W; // number of records
  int<lower=0> O_t[W]; // number of reported cases at time t
  real<lower=0> Ostar_t[W]; // exposure at time t
  int<lower=0> sumO_t[W]; // cumulative number of reported cases at time t
  int<lower=0> pop; // island population
}
parameters {
  real<lower=0> beta; 
  real<lower=0,upper=1> rho;
  real<lower=0> phi;
}
transformed parameters {
  real<lower=0> lp[W];
  real<lower=0> sampledisp[W];
  for(i in 1:W) {
    lp[i] = ( 1 - sumO_t[i] / (rho * pop)) * beta * Ostar_t[i] ;
    sampledisp[i] = lp[i]/phi;
  }
}
model {
  beta ~ exponential(0.1);
  rho ~ beta(1,1);
  phi ~ cauchy(0,2.5);
  // likelihood
  if(inference==1) target += neg_binomial_2_lpmf(O_t|lp,sampledisp);
}


