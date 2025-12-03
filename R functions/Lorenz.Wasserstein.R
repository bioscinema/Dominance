#' Wasserstein-1 distance between dominance Lorenz curves
#'
#' This function computes pairwise Wasserstein-1 distances (a.k.a. Earth
#' Mover's Distance) between samples based on their dominance Lorenz
#' curves. Given a Lorenz matrix of size \eqn{K \times n}, where rows
#' correspond to ranked taxa and columns correspond to samples, the
#' Wasserstein distance between samples \eqn{i} and \eqn{j} is:
#' \deqn{
#'   W_1(L_i, L_j) = \frac{1}{K} \sum_{k=1}^K | L_i(k) - L_j(k) |,
#' }
#' which represents the discrete integral of the absolute difference
#' between the two Lorenz curves.
#'
#' @param lorenz_matrix A numeric matrix of dimension \eqn{K \times n},
#'   where each column is a dominance Lorenz curve (e.g., output from
#'   \code{\link{dominance_lorenz}}).
#'
#' @return A \code{dist} object containing the pairwise Wasserstein-1
#'   distances between samples.
#'
#' @examples
#' set.seed(1)
#' lor <- matrix(runif(100), nrow = 20)
#' colnames(lor) <- paste0("Sample_", 1:5)
#' lor <- apply(lor, 2, function(x) cumsum(sort(x, decreasing = TRUE)) / sum(x))
#'
#' dist_w1 <- Lorenz_Wasserstein(lor)
#' dist_w1
#'
#' @export
Lorenz_Wasserstein <- function(lorenz_matrix) {
  lorenz_matrix <- as.matrix(lorenz_matrix)
  
  sample_names <- colnames(lorenz_matrix)
  n_samples    <- ncol(lorenz_matrix)
  K            <- nrow(lorenz_matrix)
  
  dist_mat <- matrix(0, nrow = n_samples, ncol = n_samples)
  colnames(dist_mat) <- rownames(dist_mat) <- sample_names
  
  for (i in 1:(n_samples - 1)) {
    for (j in (i + 1):n_samples) {
      L_i <- lorenz_matrix[, i]
      L_j <- lorenz_matrix[, j]
      
      # discrete W1 distance between Lorenz curves
      w_dist <- sum(abs(L_i - L_j)) / K
      
      dist_mat[i, j] <- dist_mat[j, i] <- w_dist
    }
  }
  
  as.dist(dist_mat)
}
