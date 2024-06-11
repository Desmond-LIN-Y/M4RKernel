data {
  int<lower=1> N;      // number of observations
  array[N] vector[2] xn ;         // univariate covariate
  vector[N] y;         // target variable
}
transformed data {
  // Normalize data
 //real xmean = mean(x);
  real ymean = mean(y);
 // real xsd = sd(x);
  real ysd = sd(y);
 // array[N] real xn = to_array_1d((x - xmean)/xsd);
  // vector[N] x1 = (x-xmean)/xsd;
  vector[N] yn = (y - ymean)/ysd;
  //array[N] real x2 = to_array_1d(ra);
  // array[2] real xn = to_array_2d(append_col(x1, ra));
  real sigma_intercept = 0.1;
  vector[N] jitter = rep_vector(1e-9, N);
}
parameters {
  vector<lower=0>[2] lengthscale_f; // lengthscale of f
  real<lower=0> sigma_f;       // scale of f
  vector<lower=0>[2] lengthscale_g; // lengthscale of g
  real<lower=0> sigma_g;       // scale of g
  vector[N] z_f;
  vector[N] z_g;
}
model {
  vector[2] mu = rep_vector(1,2);
  matrix[2,2] s2 = [[1,0],[0,1]];
  // covariances and Cholesky decompositions
  matrix[N, N] K_f = gp_exp_quad_cov(xn, sigma_f, to_array_1d(lengthscale_f))+
                     sigma_intercept;
  matrix[N, N] L_f = cholesky_decompose(add_diag(K_f, jitter));
  matrix[N, N] K_g = gp_exp_quad_cov(xn, sigma_g, to_array_1d(lengthscale_g))+
                     sigma_intercept;
  matrix[N, N] L_g = cholesky_decompose(add_diag(K_g, jitter));
  // priors
  z_f ~ std_normal();
  z_g ~ std_normal();

  lengthscale_f ~ multi_normal(mu, s2);
  lengthscale_g ~ multi_normal(mu, s2);
  sigma_f ~ normal(0, .5);
  sigma_g ~ normal(0, .5);
  // model
  yn ~ normal(L_f * z_f, exp(L_g * z_g));
}
generated quantities {
  vector[N] f;
  vector[N] sigma;
  {
    // covariances and Cholesky decompositions
    matrix[N, N] K_f = gp_exp_quad_cov(xn, sigma_f, to_array_1d(lengthscale_f))+
                       sigma_intercept;
    matrix[N, N] L_f = cholesky_decompose(add_diag(K_f, jitter));
    matrix[N, N] K_g = gp_exp_quad_cov(xn, sigma_g, to_array_1d(lengthscale_g))+
                       sigma_intercept;
    matrix[N, N] L_g = cholesky_decompose(add_diag(K_g, jitter));
    // function scaled back to the original scale
    f = (L_f * z_f)*ysd + ymean;
    sigma = exp(L_g * z_g)*ysd;
  }
}