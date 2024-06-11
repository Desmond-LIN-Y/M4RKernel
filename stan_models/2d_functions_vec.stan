vector lambda(real L, vector m) {
  return m * pi() / (2 * L); // sqrt version, m is positive
}

matrix sqrt_LAMBDA_2D(real L1, real L2, matrix S) {
  matrix[rows(S),2] sqrt_LAMBDA;
  sqrt_LAMBDA[:,1] = lambda(L1, S[,1]); // already sqrt
  sqrt_LAMBDA[:,2] = lambda(L1, S[,2]);

  return sqrt_LAMBDA;
}

matrix PHI_2D(real L1, real L2,  matrix X) {
    int N = rows(X);
    int M1 = 30;
    int M2 = 30;
    matrix[N,M1*M2] PHI;
    matrix[N,M1] PHI_1 = sin(diag_post_multiply(rep_matrix(pi()/(2*L1) * (X[,1]+L1), M1), linspaced_vector(M1, 1, M1)))/sqrt(L1);
    matrix[N,M2] PHI_2 = sin(diag_post_multiply(rep_matrix(pi()/(2*L2) * (X[,2]+L2), M2), linspaced_vector(M2, 1, M2)))/sqrt(L2);
    PHI[,1:M2] = rep_matrix(PHI_1[,1], M2) .* PHI_2;
    for(i in 2:M1)
      PHI[,1:(M2*i)] = append_col(PHI[,1:(M2*(i-1))], rep_matrix(PHI_1[,i], M2) .* PHI_2);
    return PHI;
  }



// ========= Spectral density functions =========

vector diagSPD_eq_2d(real rho1, real rho2,  real L1, real L2, matrix indices) {
    vector[2] L = to_vector([L1,L2]);
    vector[2] lscale = to_vector([rho1, rho2]);
    return  sqrt(sqrt(2*pi())^2 * rho1 * rho2) * exp(-.25 * (indices * (lscale*pi() ./ (2*L))^2));
  }


 vector diagSPD_m52_2d(real rho1, real rho2, real L1, real L2, matrix indices) {
    int J = rows(indices);
    vector[2] L = to_vector([L1,L2]);
    vector[2] lscale = to_vector([rho1, rho2]);
    return  pow(5 + (indices * (lscale*pi() ./ (2*L))^2), -7.0/4.0) * sqrt(10*pi()*rho1*rho2)*pow(5, 5.0/4.0);
  }


vector diagSPD_m32_2d(real rho1, real rho2, real L1, real L2, matrix indices) {
    int J = rows(indices);
    vector[2] L = to_vector([L1,L2]);
    vector[2] lscale = to_vector([rho1, rho2]);
  return pow(3 + (indices * (lscale*pi() ./ (2*L))^2), -5.0/4.0) * sqrt(6*pi()*rho1*rho2)*pow(3, 3.0/4.0);
}



vector diagSPD_m12_2d(real rho1, real rho2, real L1, real L2, matrix indices) {
    int J = rows(indices);
    vector[2] L = to_vector([L1,L2]);
    vector[2] lscale = to_vector([rho1, rho2]);
  return pow(1 + (indices * (lscale*pi() ./ (2*L))^2), -3.0/4.0) * sqrt(2*pi()*rho1*rho2);
}


vector hsgp_m32_2d(vector beta, real alpha, real rho1, real rho2, real L1, real L2, matrix S2, matrix PHI) {
  vector[rows(S2)] diagSPD = alpha * diagSPD_m32_2d(rho1, rho2, L1, L2, S2);
  return  PHI * (diagSPD .* beta);
}



vector hsgp_eq_2d(vector beta, real alpha, real rho1, real rho2, real L1, real L2, matrix S2, matrix PHI) {
  vector[rows(S2)] diagSPD = alpha * diagSPD_eq_2d(rho1, rho2, L1, L2, S2);
  return PHI * (diagSPD .* beta);

}

vector hsgp_m52_2d(vector beta, real alpha, real rho1, real rho2, real L1, real L2, matrix S2, matrix PHI) {
  vector[rows(S2)] diagSPD = alpha * diagSPD_m52_2d(rho1, rho2, L1, L2, S2);

  return PHI * (diagSPD .* beta);
}


vector hsgp_m12_2d(vector beta, real alpha, real rho1, real rho2, real L1, real L2,matrix S2, matrix PHI) {
  vector[rows(S2)] diagSPD = alpha * diagSPD_m12_2d(rho1, rho2, L1, L2, S2);

  return PHI * (diagSPD .* beta);
}

vector hsgp_m12m32_2d(vector beta, real alpha1, real alpha2, real rho1, real rho2, real L1, real L2,matrix S2, matrix PHI) {
  vector[rows(S2)] diagSPD = alpha1 * diagSPD_m12_2d(rho1, rho2, L1, L2, S2) + alpha2 * diagSPD_m32_2d(rho1, rho2, L1, L2, S2);

  return PHI * (diagSPD .* beta);
}
vector hsgp_m12eq_2d(vector beta, real alpha1, real alpha2, real rho1, real rho2, real L1, real L2,matrix S2, matrix PHI) {
  vector[rows(S2)] diagSPD = alpha1 * diagSPD_m12_2d(rho1, rho2, L1, L2, S2) + alpha2 * diagSPD_eq_2d(rho1, rho2, L1, L2, S2); 

  return PHI * (diagSPD .* beta);
}
vector hsgp_m32eq_2d(vector beta, real alpha1, real alpha2, real rho1, real rho2, real L1, real L2,matrix S2, matrix PHI) {
  vector[rows(S2)] diagSPD = alpha1 * diagSPD_eq_2d(rho1, rho2, L1, L2, S2) +  alpha2 * diagSPD_m32_2d(rho1, rho2, L1, L2, S2);

  return PHI * (diagSPD .* beta);
}
vector hsgp_m32m52_2d(vector beta, real alpha1, real alpha2, real rho1, real rho2, real L1, real L2,matrix S2, matrix PHI) {
  vector[rows(S2)] diagSPD = alpha1 * diagSPD_m32_2d(rho1, rho2, L1, L2, S2) +  alpha2 * diagSPD_m52_2d(rho1, rho2, L1, L2, S2);

  return PHI * (diagSPD .* beta);
}
vector hsgp_m52eq_2d(vector beta, real alpha1, real alpha2, real rho1, real rho2, real L1, real L2,matrix S2, matrix PHI) {
  vector[rows(S2)] diagSPD = alpha1 * diagSPD_m52_2d(rho1, rho2, L1, L2, S2) +  alpha2 * diagSPD_m52_2d(rho1, rho2, L1, L2, S2);

  return PHI * (diagSPD .* beta);
}vector hsgp_m52m12_2d(vector beta, real alpha1, real alpha2, real rho1, real rho2, real L1, real L2,matrix S2, matrix PHI) {
  vector[rows(S2)] diagSPD = alpha1 * diagSPD_m12_2d(rho1, rho2, L1, L2, S2) +  alpha2 * diagSPD_m52_2d(rho1, rho2, L1, L2, S2);

  return PHI * (diagSPD .* beta);
}

