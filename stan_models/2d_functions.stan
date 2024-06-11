// ========== Basis functions ==========

vector lambda(real L, vector m) {
  return m * pi() / (2 * L); // sqrt version, m is positive
}

matrix sqrt_LAMBDA_2D(real L1, real L2, matrix S) {
  matrix[rows(S),2] sqrt_LAMBDA;
  sqrt_LAMBDA[:,1] = lambda(L1, S[,1]); // already sqrt
  sqrt_LAMBDA[:,2] = lambda(L1, S[,2]);

  return sqrt_LAMBDA;
}


matrix PHI_2D(real L1, real L2, matrix S, matrix X) {
	int N = rows(X);
	int J = rows(S);

	// Precompute to avoid repetitive calculations in loop
	real inv_sqrt_L1 = 1/sqrt(L1);
	real inv_sqrt_L2 = 1/sqrt(L2);
	matrix[J,2] sqrt_LAMBDA = sqrt_LAMBDA_2D(L1, L2, S);
	vector[N] x1_plus_L1 = X[,1] + L1;
	vector[N] x2_plus_L2 = X[,2] + L2;

	matrix[N,J] PHI;
  for (j in 1:J) {
		vector[N] phi_j1 = inv_sqrt_L1 * sin(sqrt_LAMBDA[j,1] * x1_plus_L1);
		vector[N] phi_j2 = inv_sqrt_L2 * sin(sqrt_LAMBDA[j,2] * x2_plus_L2);
		PHI[,j] = phi_j1 .* phi_j2;
  }

  return PHI;
}




// Squared exponential kernel
real spd_se_2d(real alpha, real rho1, real rho2, row_vector sq_omega) {
  return alpha * 2 * pi() * rho1 * rho2 * square(exp(-0.5 * dot_product(square([rho1, rho2]'), to_vector(sq_omega))));
}

vector diagSPD_EQ_2D(real alpha, real rho1, real rho2, matrix sqrt_LAMBDA) {
  int J = rows(sqrt_LAMBDA);
  vector[J] diagSPD;
  for (j in 1:J) {
    diagSPD[j] = sqrt(spd_se_2d(alpha, rho1, rho2, square(sqrt_LAMBDA[j])));
  }

  return diagSPD;
}


// Matern 5/2 kernel
real spd_m52_2d(real alpha, real rho1, real rho2, row_vector sq_omega) {
  real d = dot_product(square([rho1, rho2]'), to_vector(sq_omega));
  return alpha * 10 * pi() * 5.0^(5.0/2.0) * rho1 * rho2 * (5 + d)^(-7.0/2.0);
}

vector diagSPD_m52_2d(real alpha, real rho1, real rho2, matrix sqrt_LAMBDA) {
  int J = rows(sqrt_LAMBDA);
  vector[J] diagSPD;
  for (j in 1:J) {
    diagSPD[j] = sqrt(spd_m52_2d(alpha, rho1, rho2, square(sqrt_LAMBDA[j])));
  }

  return diagSPD;
}

// Matern 3/2 kernel 
real spd_m32_2d(real alpha, real rho1, real rho2, row_vector sq_omega){
  real d = dot_product(square([rho1,rho2]'), to_vector(sq_omega));
  return alpha * 6 * pi() * 3.0^(3.0/2.0) * rho1 * rho2 * (3 + d)^(-5.0/2.0); 

}

vector diagSPD_m32_2d(real alpha, real rho1, real rho2, matrix sqrt_LAMBDA) {
  int J = rows(sqrt_LAMBDA);
  vector[J] diagSPD;
  for (j in 1:J) {
    diagSPD[j] = sqrt(spd_m32_2d(alpha, rho1, rho2, square(sqrt_LAMBDA[j])));
  }

  return diagSPD;
}

// Matern 1/2 kernel 
real spd_m12_2d(real alpha, real rho1, real rho2, row_vector sq_omega){
  real d = dot_product(square([rho1,rho2]'), to_vector(sq_omega));
  return alpha * 2 * pi() *  rho1 * rho2 * (1 + d)^(-3.0/2.0); 

}

vector diagSPD_m12_2d(real alpha, real rho1, real rho2, matrix sqrt_LAMBDA) {
  int J = rows(sqrt_LAMBDA);
  vector[J] diagSPD;
  for (j in 1:J) {
    diagSPD[j] = sqrt(spd_m12_2d(alpha, rho1, rho2, square(sqrt_LAMBDA[j])));
  }

  return diagSPD;
}


vector hsgp_eq_2d(vector beta, real alpha, real rho1, real rho2, matrix sqrt_LAMBDA, matrix PHI) {
  vector[rows(sqrt_LAMBDA)] diagSPD = diagSPD_EQ_2D(alpha, rho1, rho2, sqrt_LAMBDA);

  return PHI * (diagSPD .* beta);
}

vector hsgp_m52_2d(vector beta, real alpha, real rho1, real rho2, matrix sqrt_LAMBDA, matrix PHI) {
  vector[rows(sqrt_LAMBDA)] diagSPD = diagSPD_m52_2d(alpha, rho1, rho2, sqrt_LAMBDA);

  return PHI * (diagSPD .* beta);
}

vector hsgp_m32_2d(vector beta, real alpha, real rho1, real rho2, matrix sqrt_LAMBDA, matrix PHI) {
  vector[rows(sqrt_LAMBDA)] diagSPD = diagSPD_m32_2d(alpha, rho1, rho2, sqrt_LAMBDA);

  return PHI * (diagSPD .* beta);
}
vector hsgp_m12_2d(vector beta, real alpha, real rho1, real rho2, matrix sqrt_LAMBDA, matrix PHI) {
  vector[rows(sqrt_LAMBDA)] diagSPD = diagSPD_m12_2d(alpha, rho1, rho2, sqrt_LAMBDA);

  return PHI * (diagSPD .* beta);
}
