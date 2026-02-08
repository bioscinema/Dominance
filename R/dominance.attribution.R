#' Group-wise dominance summary: stable vs opportunistic dominance

#' For each taxon within each group, this function summarizes dominance
#' patterns based on a rank matrix. It computes:
#' \itemize{
#'   \item \code{Borda_norm}: normalized Borda score in [0, 1], measuring
#'         average dominance level across samples in the group;
#'   \item \code{RankEntropy_norm}: normalized rank entropy in [0, 1],
#'         measuring instability/variability of dominance ranks;
#'   \item \code{StableDominance}: \code{Borda_norm * (1 - RankEntropy_norm)},
#'         capturing taxa that are both dominant and stably ranked;
#'   \item \code{OpportunisticDominance}: \code{Borda_norm * RankEntropy_norm},
#'         capturing taxa that are dominant but fluctuate strongly across samples.
#' }
#'
#' Optionally, taxonomy information can be merged if a \code{phyloseq}
#' object is provided.
#'
#' @param contrib_matrix A numeric matrix of taxa x samples containing contrib
#'   information (e.g., output from \code{\link{Contrib.eval}}), where
#'   smaller ranks correspond to more dominant taxa. Rows are taxa and
#'   columns are samples.
#' @param metadata A \code{data.frame} with row names matching the sample
#'   IDs (column names of \code{rank_matrix}), containing grouping
#'   variables.
#' @param group_var Character string; the column name in \code{metadata}
#'   specifying the grouping factor.
#' @param physeq Optional \code{phyloseq} object. If provided, the taxonomy
#'   table (\code{tax_table(physeq)}) will be joined to the output by taxon
#'   name.
#'
#' @details
#' For a given group \eqn{g} with \eqn{n_g} samples and \eqn{K} taxa:
#' \itemize{
#'   \item \strong{Normalized Borda}:
#'     \deqn{
#'       \text{Borda}_i = \sum_{s \in g} (K - \text{rank}_{i,s}), \quad
#'       \text{Borda\_norm}_i =
#'         \frac{\text{Borda}_i}{n_g (K - 1)} \in [0, 1],
#'     }
#'     where values close to 1 indicate consistently high dominance ranks.
#'
#'   \item \strong{Normalized rank entropy}:
#'     Let \eqn{p_{i,r}} be the empirical frequency of rank \eqn{r} for taxon
#'     \eqn{i} within group \eqn{g}, and
#'     \deqn{
#'       H_i = -\sum_r p_{i,r}\log_2(p_{i,r})
#'     }
#'     be the Shannon entropy of the rank distribution. The maximum possible
#'     entropy is approximated as:
#'     \deqn{
#'       H_{\max} = \log_2(\min(K, n_g)),
#'     }
#'     and the normalized rank entropy is
#'     \deqn{
#'       \text{RankEntropy\_norm}_i =
#'         \begin{cases}
#'           H_i / H_{\max}, & H_{\max} > 0, \\
#'           0, & H_{\max} = 0.
#'         \end{cases}
#'     }
#'     Values near 0 indicate stable ranks across samples, while values near
#'     1 indicate highly variable ranks.
#'
#'   \item \strong{Stable dominance}:
#'     \deqn{
#'       \text{StableDominance}_i = \text{Borda\_norm}_i
#'                                 \cdot (1 - \text{RankEntropy\_norm}_i),
#'     }
#'     high when a taxon is both strongly and stably dominant.
#'
#'   \item \strong{Opportunistic dominance}:
#'     \deqn{
#'       \text{OpportunisticDominance}_i = \text{Borda\_norm}_i
#'                                        \cdot \text{RankEntropy\_norm}_i,
#'     }
#'     high when a taxon is often dominant but with highly fluctuating ranks,
#'     resembling opportunistic or bloom-like behaviour.
#' }
#'
#' @return A \code{data.frame} with one row per taxon per group, containing:
#'   \itemize{
#'     \item \code{Taxon}: taxon ID (row name of \code{rank_matrix});
#'     \item \code{Group}: group label;
#'     \item \code{Borda_norm}: normalized Borda score in [0, 1];
#'     \item \code{RankEntropy_norm}: normalized rank entropy in [0, 1];
#'     \item \code{StableDominance}: stable dominance index;
#'     \item \code{OpportunisticDominance}: opportunistic dominance index;
#'     \item additional taxonomy columns if \code{physeq} is provided.
#'   }
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#'
#' set.seed(1)
#' K <- 10; n <- 20
#' rank_mat <- matrix(
#'   sample(1:K, K * n, replace = TRUE),
#'   nrow = K, ncol = n
#' )
#' rownames(rank_mat) <- paste0("Taxon_", 1:K)
#' colnames(rank_mat) <- paste0("Sample_", 1:n)
#'
#' meta <- data.frame(
#'   Group = rep(c("A", "B"), each = 10),
#'   row.names = colnames(rank_mat)
#' )
#'
#' res <- EfficientDominanceSummary(
#'   rank_matrix = rank_mat,
#'   metadata = meta,
#'   group_var = "Group"
#' )
#'
#' head(res)
#' }
#' @importFrom phyloseq tax_table
#' @importFrom magrittr %>%
#' @export
EfficientDominanceSummary <- function(contrib_matrix, metadata, group_var, physeq = NULL) {

  Dominance.rank <- function(contrib_matrix) {

    rank_matrix <- apply(contrib_matrix, 2, function(x) {
      rank(-x, ties.method = "average")
    })

    dimnames(rank_matrix) <- dimnames(contrib_matrix)

    return(rank_matrix)
  }

  rank_matrix <- Dominance.rank(contrib_matrix)
  rank_matrix <- as.matrix(rank_matrix)

  # Align metadata to columns (samples) of rank_matrix
  metadata <- metadata[match(colnames(rank_matrix), rownames(metadata)), ]
  stopifnot(all(rownames(metadata) == colnames(rank_matrix)))

  groups      <- unique(metadata[[group_var]])
  K           <- nrow(rank_matrix)
  taxon_names <- rownames(rank_matrix)

  taxonomy_df <- NULL
  if (!is.null(physeq)) {
    taxonomy_df <- as.data.frame(tax_table(physeq))
    taxonomy_df$Taxon <- rownames(taxonomy_df)
  }

  result_list <- lapply(groups, function(g) {
    samples_in_group <- rownames(metadata)[metadata[[group_var]] == g]
    R_sub <- rank_matrix[, samples_in_group, drop = FALSE]

    n_g <- ncol(R_sub)

    # ----- Normalized Borda -----
    # Borda_raw_i = sum_s (K - rank_{i,s}), Borda_norm in [0,1]
    borda_raw  <- rowSums(K - R_sub)
    if (K > 1 && n_g > 0) {
      borda_norm <- borda_raw / (n_g * (K - 1))
    } else {
      borda_norm <- rep(NA_real_, K)
    }

    # ----- Normalized rank entropy -----
    # For each taxon, entropy of rank distribution over samples
    shannon <- apply(R_sub, 1, function(ranks) {
      freq <- tabulate(ranks, nbins = K) / length(ranks)
      freq <- freq[freq > 0]
      if (length(freq) == 0) return(0)
      -sum(freq * log2(freq))
    })

    H_max <- if (n_g > 0) log2(min(K, n_g)) else 0
    if (H_max > 0) {
      RankEntropy_norm <- shannon / H_max
    } else {
      RankEntropy_norm <- rep(0, K)
    }

    # Truncate to [0,1] for numerical safety
    RankEntropy_norm <- pmin(pmax(RankEntropy_norm, 0), 1)
    borda_norm       <- pmin(pmax(borda_norm,       0), 1)

    # ----- Composite indices -----
    StableDominance        <- borda_norm * (1 - RankEntropy_norm)
    OpportunisticDominance <- borda_norm * RankEntropy_norm

    df <- data.frame(
      Taxon                 = taxon_names,
      Group                 = g,
      Borda_norm            = borda_norm,
      RankEntropy_norm      = RankEntropy_norm,
      StableDominance       = StableDominance,
      OpportunisticDominance = OpportunisticDominance,
      stringsAsFactors      = FALSE
    )

    if (!is.null(taxonomy_df)) {
      df <- dplyr::left_join(df, taxonomy_df, by = "Taxon")
    }

    df %>% dplyr::arrange(dplyr::desc(StableDominance))
  })

  result <- dplyr::bind_rows(result_list)
  result
}
