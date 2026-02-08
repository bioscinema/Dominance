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
#' Compute dominance Lorenz matrix from an OTU table
#'
#' @param otu A numeric matrix of counts or abundances
#'   (taxa in rows, samples in columns).
#' @param method Character string specifying the contribution kernel;
#'   one of "Shannon", "Simpson", or "identity".
#' @param TSS Logical; if TRUE (default), each sample is scaled
#'   to sum to 1 (total-sum scaling).
#' @return A numeric matrix (taxa x samples) of dominance Lorenz values.
#' @export
dominance_lorenz <- function(otu,
                             method = c("Shannon", "Simpson", "identity"),
                             TSS = TRUE) {
  method <- match.arg(method)
  otu    <- as.matrix(otu)

  Lorenz.eval <- function(contrib_matrix) {
    # input: taxon x sample
    # output: taxon x sample

    lorenz_matrix <- matrix(NA, nrow = nrow(contrib_matrix), ncol = ncol(contrib_matrix))
    colnames(lorenz_matrix) <- colnames(contrib_matrix)

    for (j in seq_len(ncol(contrib_matrix))) {
      contrib_j <- contrib_matrix[, j]
      contrib_sorted <- sort(contrib_j, decreasing = TRUE)
      cum_contrib <- cumsum(contrib_sorted)
      total_contrib <- sum(contrib_sorted)
      lorenz_matrix[, j] <- cum_contrib / total_contrib
    }

    return(lorenz_matrix)
  }

  if (is.null(colnames(otu))) {
    colnames(otu) <- paste0("Sample_", seq_len(ncol(otu)))
  }
  if (is.null(rownames(otu))) {
    rownames(otu) <- paste0("Taxon_", seq_len(nrow(otu)))
  }

  rel_abun <- otu

  if (TSS) {
    col_sum  <- colSums(otu, na.rm = TRUE)
    zero_col <- (col_sum == 0)
    non_zero <- !zero_col

    if (any(zero_col)) {
      warning("Some samples have total count 0; their relative abundances are set to 0.")
    }

    if (any(non_zero)) {
      rel_abun[, non_zero] <- sweep(
        otu[, non_zero, drop = FALSE],
        2,
        col_sum[non_zero],
        "/"
      )
    }
    if (any(zero_col)) {
      rel_abun[, zero_col] <- 0
    }
  }

  contrib <- Contrib.eval(rel_abun, method = method)
  lorenz  <- Lorenz.eval(contrib)

  colnames(lorenz) <- colnames(otu)
  lorenz
}
