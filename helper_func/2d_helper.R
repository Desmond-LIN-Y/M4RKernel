make_S_matrix <- function(M1, M2){# this creates the indexing matrix required in spectral density function
  S <- expand.grid(seq(M1), seq(M2))
  S <- as.matrix(S)
  tmp <- S[,2]
  S[,2] <- S[,1]
  S[,1] <- tmp
  
  return(S)
}

make_lowertri_idxset <- function(A) {
  x <- matrix(1:A*A, nrow = A, ncol = A)
  
  # Indices of the lower triangular part, excluding the diagonal
  lower_tri_indices <- which(row(x) >= col(x), arr.ind = TRUE)
  
  return(lower_tri_indices)
}
make_nn_idxset <- function(A) {
  nn_indices <- matrix(NA, nrow = A*A, ncol = 2)
  n <- 1
  for (i in 1:A) {
    for (j in 1:A) {
      nn_indices[n,1] <- A - i + j
      nn_indices[n,2] <- i
      
      n <- n + 1
    }
  }
  
  return(nn_indices)
}


make_nn_lowertri_idxset <- function(A, parametrization='rd'){
  # Make an index set for a full square matrix
  dt_full_idxset <- expand.grid(key1 = 1:A, key2 = 1:A)
  
  # Make an non-nuisance index set
  if (parametrization == 'rd'){
  dt_nn_idxset <- as.data.frame(make_nn_idxset(A))}
  else{dt_nn_idxset <- as.data.frame(list(
    V1 = rep(1:A, times=A),
    V2 = rep(1:A, each=A)
  ), col.names = c("V1","V2"))}
  
  # Make a lower triangular index set
  dt_lowertri_idxset <- as.data.frame(make_lowertri_idxset(A))
  
  # Combine and sort to produce lower triangular index set for a nn index set
  colnames(dt_lowertri_idxset) <- c("key1", "key2")
  dt_nn_lowertri_idxset <- merge(dt_lowertri_idxset, cbind(dt_full_idxset, dt_nn_idxset), by = c("key1", "key2"))
  dt_nn_lowertri_idxset <- dplyr::arrange(dt_nn_lowertri_idxset, V2, V1)
  dt_nn_lowertri_idxset <- dplyr::select(dt_nn_lowertri_idxset, V1, V2)
  
  return(as.matrix(dt_nn_lowertri_idxset))
}
make_sym_from_lowertri_idx <- function(A) { # could be used for rectangular matrices
  index_set <- rep(NA, A)
  n <- 1
  
  for (j in 1:A) {
    for (i in 1:A) {
      if (i >= j) {
        index_set[n] <- i + (j - 1)*A - j*(j - 1)/2
      } else {
        index_set[n] <- j + (i - 1)*A - i*(i - 1)/2
      }
      
      n <- n + 1
    }
  }
  
  return(index_set)
}