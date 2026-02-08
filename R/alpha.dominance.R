#' Alpha-scale dominance decomposition: scarcity + Gini
#'
#' Decompose within-sample dominance into:
#' \itemize{
#'   \item \code{scarcity}: proportion of non-dominant taxa (with optional Chao-style correction)
#'   \item \code{Gini}: Gini index among dominant taxa
#' }
#'
#' @param otu A numeric matrix of taxa x samples (count or relative abundance).
#' @param method Contribution kernel passed to \code{\link{Contrib.eval}}:
#'   one of \code{"Shannon"}, \code{"Simpson"}, \code{"identity"}.
#' @param dominance_threshold Numeric in (0,1], cumulative contribution threshold for dominance set.
#' @param abundance_type One of \code{"count"} or \code{"relative"}.
#' @param count_tolerance Tolerance for estimating f1/f2 from relative abundance.
#' @param correction Logical; if \code{TRUE} (default), apply Chao-style bias correction
#'   to raw emptiness. If \code{FALSE}, use raw emptiness only.
#'
#' @return A \code{data.frame} with rownames = sample names and columns:
#'   \itemize{
#'     \item \code{scarcity} — proportion of non-dominant taxa
#'     \item \code{Gini} — conditional dominance Gini index
#'   }
#'
#' @export
alpha_dominance <- function(otu,
                            method = c("Shannon", "Simpson", "identity"),
                            dominance_threshold = 1,
                            abundance_type = c("count", "relative"),
                            count_tolerance = 1e-6,
                            correction = TRUE) {

  method         <- match.arg(method)
  abundance_type <- match.arg(abundance_type)

  otu <- as.matrix(otu)
  if (is.null(colnames(otu))) {
    colnames(otu) <- paste0("Sample_", seq_len(ncol(otu)))
  }
  if (is.null(rownames(otu))) {
    rownames(otu) <- paste0("Taxon_", seq_len(nrow(otu)))
  }

  K        <- nrow(otu)
  samples  <- colnames(otu)

  # ----- Convert to relative abundance if needed -----
  if (abundance_type == "count") {
    col_sum  <- colSums(otu, na.rm = TRUE)
    rel_abun <- matrix(0, nrow = nrow(otu), ncol = ncol(otu),
                       dimnames = dimnames(otu))
    non_zero <- col_sum > 0
    if (any(non_zero)) {
      rel_abun[, non_zero] <- sweep(
        otu[, non_zero, drop = FALSE],
        2, col_sum[non_zero], "/"
      )
    }
  } else {
    rel_abun <- otu
    col_sum  <- colSums(rel_abun, na.rm = TRUE)
    non_zero <- col_sum > 0
    if (any(non_zero)) {
      rel_abun[, non_zero] <- sweep(
        rel_abun[, non_zero, drop = FALSE],
        2, col_sum[non_zero], "/"
      )
    }
  }

  # ----- Contribution matrix -----
  contrib_matrix <- Contrib.eval(rel_abun, method = method)

  # ----- Output container -----
  results <- data.frame(
    scarcity = rep(NA_real_, length(samples)),
    Gini     = rep(NA_real_, length(samples)),
    row.names = samples
  )

  # ----- Loop over samples -----
  for (j in seq_len(ncol(contrib_matrix))) {
    contrib_j <- contrib_matrix[, j]

    if (all(is.na(contrib_j)) || sum(contrib_j) <= 0) {
      results$scarcity[j] <- NA_real_
      results$Gini[j]     <- NA_real_
      next
    }

    contrib_sorted <- sort(contrib_j, decreasing = TRUE)
    total          <- sum(contrib_sorted)

    cum_contrib <- cumsum(contrib_sorted)
    m <- which(cum_contrib / total >= dominance_threshold)[1]

    raw_emptiness <- K - m

    # ----- Compute f1 and f2 if correction = TRUE -----
    if (correction) {
      if (abundance_type == "count") {
        counts_j <- otu[, j]
        f1 <- sum(counts_j == 1)
        f2 <- sum(counts_j == 2)
      } else {
        rel_j <- rel_abun[, j]
        positive <- rel_j[rel_j > 0]
        if (length(positive) == 0) {
          f1 <- 0; f2 <- 0
        } else {
          N_hat <- round(1 / min(positive))
          if (!is.finite(N_hat) || N_hat <= 0) {
            f1 <- 0; f2 <- 0
          } else {
            f1 <- sum(abs(rel_j - 1 / N_hat) < count_tolerance)
            f2 <- sum(abs(rel_j - 2 / N_hat) < count_tolerance)
          }
        }
      }

      if (f2 > 0) {
        correction_term <- (f1^2) / (2 * f2)
      } else if (f1 > 1) {
        correction_term <- (f1 * (f1 - 1)) / (2 * (f2 + 1))
      } else {
        correction_term <- 0
      }

      emptiness_bc <- raw_emptiness + correction_term

    } else {
      # No correction
      emptiness_bc <- raw_emptiness
    }

    # Scale to [0, 1]
    emptiness_bc <- max(0, min(emptiness_bc, K))
    scarcity     <- emptiness_bc / K

    # ----- Conditional Gini -----
    if (m <= 1) {
      cond_gini <- 0
    } else {
      p_top <- contrib_sorted[1:m]
      if (sum(p_top) <= 0) {
        cond_gini <- 0
      } else {
        p_top_norm <- p_top / sum(p_top)
        cum_p      <- cumsum(p_top_norm)
        t_seq      <- seq_len(m) / m
        cond_gini  <- sum(cum_p - t_seq) * (2 / m)
      }
    }

    results$scarcity_ratio[j] <- scarcity
    results$Gini[j]     <- cond_gini
    results$scarcity[j] <- emptiness_bc
  }

  results
}
