functions
{
real spd_2D(real alpha, real rho1, real rho2, real w1, real w2) {
		real S;
		S = alpha^2 * sqrt(2*pi())^2 * rho1*rho2 * exp(-0.5*(rho1^2*w1^2 + rho2^2*w2^2));
				
		return S;
	}

 matrix PHI_2D(int N, int M1, int M2, real L1, real L2, vector x1, vector x2) {
    matrix[N,M1*M2] PHI;
    matrix[N,M1] PHI_1 = sin(diag_post_multiply(rep_matrix(pi()/(2*L1) * (x1+L1), M1), linspaced_vector(M1, 1, M1)))/sqrt(L1);
    matrix[N,M2] PHI_2 = sin(diag_post_multiply(rep_matrix(pi()/(2*L2) * (x2+L2), M2), linspaced_vector(M2, 1, M2)))/sqrt(L2);
    PHI[,1:M2] = rep_matrix(PHI_1[,1], M2) .* PHI_2;
    for(i in 2:M1)
      PHI[,1:(M2*i)] = append_col(PHI[,1:(M2*(i-1))], rep_matrix(PHI_1[,i], M2) .* PHI_2);
    return PHI;
  }

  vector sqrt_diagSPD_2D(real gpscale, vector lscale, vector L, matrix indices, int D) {
    return gpscale *  sqrt(sqrt(2*pi())^D * prod(lscale)) * exp(-.25 * (indices.^2 * (lscale*pi() ./ (2*L))^2));
  }
matrix hsgp_spd_eq(int A, real alpha, real rho1, real rho2, real L1, real L2, int M1, int M2, matrix indices, matrix PHI, vector z){
vector[M1*M2] SPD = z .* sqrt_diagSPD_2D(alpha, to_vector([rho1,rho2]), to_vector([L1,L2]), indices, 2);
matrix[A, M1*M2] f = diag_post_multiply(PHI, SPD);
matrix[A, A] f_restruct = f[,1:A];
return f_restruct;
}

}
data{
int<lower=10> A;
 vector[A] age_idx_std;         // Standardized age index
  real<lower=0> C1; // Factor to determine the boundary value L (cohort age dimension)
  int<lower=1> M1;  // Number of basis functions (cohort age dimension)
  real<lower=0> C2; // Factor to determine the boundary value L for age of contacted individuals (age difference dimension)
  int<lower=1> M2;  // Number of basis functions (age difference dimension)
  matrix[M1*M2, 2] indices;
}
transformed data{
  real L1, L2;
  matrix[A,M1*M2] PHI;
  L1 = C1 * max(age_idx_std);
  L2 = C2 * max(age_idx_std);
  PHI = PHI_2D(A, M1,M2,L1,L2 ,age_idx_std,age_idx_std);
  real alpha = 0.2;
  real rho1 = 1;
  real rho2 = 1;
}

generated quantities {
  vector[M1*M2] z;
  matrix[A, A] f;

  {
    z = multi_normal_cholesky_rng(rep_vector(0, M1*M2), diag_matrix(rep_vector(1, M1*M2)));
    f = hsgp_spd_eq(A, alpha, rho1, rho2, L1, L2, M1, M2, indices, PHI, z);
    f = symmetrize_from_lower_tri(f);
  }
}
