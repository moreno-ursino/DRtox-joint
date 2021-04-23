data {
  int n_pat;
  int nb_doses;
  int y[n_pat];
  int dose_donnees[n_pat];
  vector[nb_doses] wm;
}
transformed data {
  vector[nb_doses] doses_trans;
  doses_trans = log(wm ./ (1-wm));
}
parameters {
  real beta0;
  real<lower=0> beta1;
}
transformed parameters { 
  vector[nb_doses] ptox;
  ptox = exp(beta0+beta1*doses_trans) ./ (1+exp(beta0+beta1*doses_trans));
}
model {
  beta0 ~ normal(0, sqrt(10));
  beta1 ~ exponential(1);
  for(i in 1:n_pat){
    y[i] ~ bernoulli(ptox[dose_donnees[i]]);
  }
}

