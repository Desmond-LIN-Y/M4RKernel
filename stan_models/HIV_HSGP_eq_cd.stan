functions
{
#include gp-functions.stan
}
data
{
  int<lower=1> A1; // length of participants' age: 15-49 expanded to 15-69
  int<lower=1> A2; // length of contacts' age: 15-69
  int<lower=1> Nc; // number of communities

  // array[A1] real age1; // age of participants
  // array[A2] real age2; // age of contacts
 vector[A1] age_std; // standardised age of contacting individuals
// standardised age of contacted individuals


  int<lower=0> Nmf; // number of observed male-female contacts for all communities
  int<lower=0> Nfm; // number of observed female-male contacts for all communities
  int<lower=0> Nmf_all; // number of observed female-male contacts for all communities including 0 part

  array[Nmf] int<lower=1> ymf_idx; // idx of mf without 0 part
  array[Nfm] int<lower=1> yfm_idx;

  array[Nmf] int<lower=0> ymf; // the reported male-female contacts
  array[Nfm] int<lower=0> yfm; // the reported female-male contacts

  array[Nmf_all] int<lower=1> ymf_rowmajor_matrix_index_all; // the row-major index in the full contact matrix
  array[Nmf_all] int<lower=1> yfm_rowmajor_matrix_index_all; // the row-major index in the full contact matrix

  array[A1*A2] int<lower=1> index_mf; // the row-major index in the single contact matrix
  array[A1*A2] int<lower=1> index_fm; // the row-major index in the single contact matrix

  vector[Nfm] log_offset_fm; // the offset terms excluding the pop sizes
  vector[Nmf] log_offset_mf;
  vector[Nmf_all] log_pop_fm_all; // the pop terms
  vector[Nmf_all] log_pop_mf_all;

  array[Nmf_all] int<lower=1> community_index_mf_all;
  array[Nmf_all] int<lower=1> community_index_fm_all;

  // HSGP arguments
  real<lower=0> c_age1; // factor c to determine the boundary value L for age
  int<lower=1> M_age1; // number of basis functions
  real<lower=0> c_age2; // factor c to determine the boundary value L for age
  int<lower=1> M_age2; // number of basis functions



}

transformed data
{
  int N = Nmf + Nfm;
  int N_all = Nmf_all + Nmf_all;
  int<lower=1> A_squared = A2 * A2;

  array[N] int y = append_array( ymf, yfm);
  real L_age1;
  real L_age2;

  // HSGP basis functions on 2D surface

  L_age2 = c_age2*max(age_std);
  L_age1 = c_age1*max(age_std);

  // Precompute HSGP basis functions
  matrix[A2, M_age1] PHI1 = PHI(A2, M_age1, L_age1, age_std);
  matrix[A1, M_age2] PHI2 = PHI(A1, M_age2, L_age2, age_std);
  array[Nmf] int<lower=1> ymf_rowmajor_matrix_index = ymf_rowmajor_matrix_index_all[ymf_idx]; // the row-major index of the observations in the full contact matrix
  array[Nfm] int<lower=1> yfm_rowmajor_matrix_index = yfm_rowmajor_matrix_index_all[yfm_idx]; // the row-major index of the observations in the full contact matrix

  array[Nmf] int<lower=1> community_index_mf = community_index_mf_all[ymf_idx]; // index for the community of observations
  array[Nfm] int<lower=1> community_index_fm = community_index_fm_all[yfm_idx];
}

parameters
{
  vector[Nc] log_effect_community; // Baseline effects
  // update for this parameter: so that for old participants, the under-reporting effects would be converaged to 0
  // real log_effect_under_report_baseline_f; // Baseline for under-reporting effect
  real<lower=0> overdispersion;

  real<lower=0> lscale_1; // GP hyperparameters
  real<lower=0> lscale_2;
  real<lower=0> gpscale;
  matrix[M_age1, M_age2] z; // GP variables
  matrix[M_age1, M_age2] z_ur;
  real<lower=0> lscale_ur_1; // GP hyperparameters of under-reporting
  real<lower=0> lscale_ur_2;
  real<lower=0> gpscale_ur;
}

transformed parameters
{
  vector[A_squared] f_fm; // Ajusted GP approx for the shared contact rates
  vector[A_squared] log_effect_age_under_report_f; // Adjusted GP approx for female under-reported effects

  // log of the contact intensities
  vector[Nmf_all] log_m_mf;
  vector[Nmf_all] log_m_fm;

  // contact intensities including age-heaped effects --> related to observations
  vector[Nmf_all] log_under_reported_m_fm;

  //age-specific effect to account for the under-reporting for females
  {
   matrix[A1, A2] raw_f; // GP approx for f_fm
   matrix[A1, A2] raw_log_effect_ur_f; // GP approx for log_effect_age_under_report_f

    // local scope
    raw_f =  hsgp(A1, gpscale, lscale_1, lscale_2, L_age1, L_age2, M_age1, M_age2, PHI1, PHI2, z);
    raw_log_effect_ur_f =  hsgp(A2, gpscale_ur, lscale_ur_1, lscale_ur_2, L_age1, L_age2, M_age1, M_age2, PHI1, PHI2, z_ur);

    // remove the extra over-reported part of women to men (including the baseline of the women effects)
    f_fm = to_vector(raw_f);
    log_effect_age_under_report_f = to_vector(raw_log_effect_ur_f);
    

    // log of the contact intensities from the model
    log_m_mf =
      log_effect_community[community_index_mf_all] + f_fm[ymf_rowmajor_matrix_index_all] + log_pop_mf_all;
    log_m_fm =
      log_effect_community[community_index_fm_all] + f_fm[yfm_rowmajor_matrix_index_all] + log_pop_fm_all;


    // age-heaped contact intensities to fit the observations
    // women to men, including the under-reporting effects with the baseline
    log_under_reported_m_fm = log_m_fm + log_effect_age_under_report_f[yfm_rowmajor_matrix_index_all];

  }
}

model
{
  // GP priors
  //
  target += inv_gamma_lpdf(lscale_1 | 5, 5);
  target += inv_gamma_lpdf(lscale_2 | 5, 5);
  target += cauchy_lpdf(gpscale | 0, 2);
  target += std_normal_lpdf( to_vector(z));

  // 2DGP priors of females
  target += inv_gamma_lpdf(lscale_ur_1 | 5, 5);
  target += inv_gamma_lpdf(lscale_ur_2 | 5, 5);
  target += cauchy_lpdf(gpscale_ur | 0, 2);
  target += std_normal_lpdf(to_vector(z_ur));

  // overdispersion
  // 240311, update the mean as 0.1 rather than 0.01
  target += exponential_lpdf(overdispersion | 10);

  // baseline for women
  // target += normal_lpdf(log_effect_under_report_baseline_f | 0, 5);

  // baseline at community level
  target += normal_lpdf(log_effect_community | 0, 10);

 {
  // combine mf and fm
  vector[N] log_under_reported_mu =
      append_row(
        log_m_mf[ymf_idx] + log_offset_mf
        ,
        log_under_reported_m_fm[yfm_idx] + log_offset_fm
      );

  target += neg_binomial_lpmf( y | exp(log_under_reported_mu + 1e-13) / overdispersion, inv(overdispersion));
  }
}

generated quantities
{
  array[Nmf] int y_pred_mf;
  array[Nfm] int y_pred_fm;

  {
   vector[Nmf] log_mu_mf_pred = log_m_mf[ymf_idx] + log_offset_mf;
   vector[Nfm] log_mu_fm_pred = log_under_reported_m_fm[yfm_idx] + log_offset_fm;

   y_pred_mf = neg_binomial_rng( exp(log_mu_mf_pred + 1e-13) / overdispersion, inv(overdispersion));
   y_pred_fm = neg_binomial_rng( exp(log_mu_fm_pred + 1e-13 ) / overdispersion, inv(overdispersion));
  }
}
