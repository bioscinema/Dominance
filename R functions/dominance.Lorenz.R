## ----------------- Helper: Hconf -----------------
#' Piecewise entropy-like contribution function
#'
#' This function maps a relative abundance \eqn{x \in [0, 1]} to an
#' entropy-like contribution. It is symmetric around \eqn{1/e} and
#' bounded by a constant \eqn{2 \log_2(e) / e}.
#'
#' @param x A numeric scalar in [0, 1].
#'
#' @return A numeric scalar giving the contribution of \code{x}.
#' @keywords internal
Hconf <- function(x) {
  cont0 <- 2 * log2(exp(1)) / exp(1)
  
  if (x == 0) {
    return(0)
  } else if (x == 1) {
    return(cont0)
  } else {
    H <- -x * log2(x)
    if (x <= 1 / exp(1)) {
      return(H)
    } else {
      return(cont0 - H)
    }
  }
}

## ----------------- Contrib.eval -----------------
#' Evaluate taxon-wise contributions for dominance analysis
#'
#' Given a matrix of relative abundances (taxa by samples), compute
#' taxon-wise contributions under different dominance kernels.
#'
#' @param rel_abun A numeric matrix of relative abundances
#'   (taxa in rows, samples in columns).
#' @param method Character string specifying the contribution kernel.
#'   One of \code{"Shannon"}, \code{"Simpson"}, or \code{"identity"}:
#'   \itemize{
#'     \item \code{"Simpson"}: squared relative abundance, \eqn{p^2}.
#'     \item \code{"Shannon"}: entropy-like transform via \code{\link{Hconf}}.
#'     \item \code{"identity"}: use the relative abundance itself, \eqn{p}.
#'   }
#'
#' @return A numeric matrix of the same dimension as \code{rel_abun},
#'   containing the contribution of each taxon in each sample.
#' @export
Contrib.eval <- function(rel_abun,
                         method = c("Shannon", "Simpson", "identity")) {
  method  <- match.arg(method)
  rel_abun <- as.matrix(rel_abun)
  
  contrib <- switch(
    method,
    "Simpson" = {
      rel_abun^2
    },
    "Shannon" = {
      apply(rel_abun, c(1, 2), Hconf)
    },
    "identity" = {
      rel_abun
    }
  )
  
  dimnames(contrib) <- dimnames(rel_abun)
  contrib
}

## ----------------- Lorenz.eval -----------------
#' Construct dominance Lorenz matrix from contributions
#'
#' For each sample, taxa are sorted by their contribution in decreasing
#' order and cumulative sums are normalized to obtain a Lorenz-type
#' dominance curve.
#'
#' @param contrib_matrix A numeric matrix of taxon-wise contributions
#'   (taxa in rows, samples in columns).
#'
#' @return A numeric matrix of the same dimension as \code{contrib_matrix},
#'   whose columns represent the dominance Lorenz curve evaluated over
#'   the sorted taxa.
#' @export
Lorenz.eval <- function(contrib_matrix) {
  contrib_matrix <- as.matrix(contrib_matrix)
  
  lorenz_matrix <- matrix(
    NA_real_,
    nrow = nrow(contrib_matrix),
    ncol = ncol(contrib_matrix)
  )
  colnames(lorenz_matrix) <- colnames(contrib_matrix)
  rownames(lorenz_matrix) <- rownames(contrib_matrix)
  
  for (j in seq_len(ncol(contrib_matrix))) {
    contrib_j      <- contrib_matrix[, j]
    contrib_sorted <- sort(contrib_j, decreasing = TRUE)
    total_contrib  <- sum(contrib_sorted)
    
    if (total_contrib == 0) {
      # If the whole column is zero, define the Lorenz curve as all zeros.
      lorenz_matrix[, j] <- 0
    } else {
      cum_contrib <- cumsum(contrib_sorted)
      lorenz_matrix[, j] <- cum_contrib / total_contrib
    }
  }
  
  lorenz_matrix
}

## ----------------- One-stop: dominance_lorenz -----------------
#' Compute dominance Lorenz matrix from an OTU table
#'
#' This function takes an OTU table (or any taxa-by-samples count matrix),
#' optionally performs total-sum scaling (TSS) to obtain relative abundances,
#' transforms them into contributions using \code{\link{Contrib.eval}},
#' and finally constructs dominance Lorenz curves via \code{\link{Lorenz.eval}}.
#'
#' @param otu A numeric matrix of counts or abundances
#'   (taxa in rows, samples in columns).
#' @param method Character string specifying the contribution kernel;
#'   one of \code{"Shannon"}, \code{"Simpson"}, or \code{"identity"}.
#'   See \code{\link{Contrib.eval}} for details.
#' @param TSS Logical; if \code{TRUE} (default), each sample is scaled
#'   to sum to 1 (total-sum scaling) before computing contributions.
#'
#' @return A numeric matrix (taxa x samples) containing the dominance
#'   Lorenz values for each taxon rank in each sample. Columns correspond
#'   to samples; within each column, rows are ordered by decreasing
#'   contribution.
#'
#' @examples
#' set.seed(1)
#' otu <- matrix(rpois(50, lambda = 10), nrow = 10, ncol = 5)
#' rownames(otu) <- paste0("Taxon_", seq_len(10))
#' colnames(otu) <- paste0("Sample_", seq_len(5))
#'
#' # Shannon-type dominance Lorenz curves
#' L_shannon  <- dominance_lorenz(otu, method = "Shannon",  TSS = TRUE)
#'
#' # Simpson-type dominance Lorenz curves
#' L_simpson  <- dominance_lorenz(otu, method = "Simpson",  TSS = TRUE)
#'
#' # Identity kernel (pure relative abundance)
#' L_identity <- dominance_lorenz(otu, method = "identity", TSS = TRUE)
#' @export
dominance_lorenz <- function(otu,
                             method = c("Shannon", "Simpson", "identity"),
                             TSS = TRUE) {
  method <- match.arg(method)
  otu    <- as.matrix(otu)
  
  rel_abun <- otu
  
  if (TSS) {
    col_sum  <- colSums(otu, na.rm = TRUE)
    zero_col <- (col_sum == 0)
    
    if (any(zero_col)) {
      warning("Some samples have total count 0; their relative abundances are set to 0.")
    }
    
    if (any(!zero_col)) {
      rel_abun[, !zero_col, drop = FALSE] <- sweep(
        otu[, !zero_col, drop = FALSE],
        2,
        col_sum[!zero_col],
        "/"
      )
    }
    
    if (any(zero_col)) {
      rel_abun[, zero_col] <- 0
    }
  }
  
  contrib <- Contrib.eval(rel_abun, method = method)
  lorenz  <- Lorenz.eval(contrib)
  
  # Keep sample names consistent with the original OTU table
  colnames(lorenz) <- colnames(otu)
  lorenz
}
