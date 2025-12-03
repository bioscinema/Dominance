#' Polar ordination along a dominance (Gini) gradient
#'
#' This function constructs a 1D ordination axis that interpolates
#' samples between two poles defined by minimum and maximum Gini
#' dominance. Given a distance object (or symmetric distance matrix)
#' and the output of \code{\link{alpha_dominance}}, it:
#' \enumerate{
#'   \item aligns samples between the distance matrix and the
#'         alpha-dominance summary;
#'   \item identifies the samples with minimum and maximum \code{Gini}
#'         as the two poles (axis=0 and axis=1);
#'   \item maps each sample to a position between 0 and 1 based on
#'         its distance to the minimum-Gini pole, normalized by the
#'         distance between the two poles.
#' }
#'
#' @param dist_obj A \code{dist} object or a symmetric numeric matrix
#'   containing pairwise distances between samples (e.g., output from
#'   \code{\link{Lorenz_Wasserstein}}).
#' @param alpha_dom A \code{data.frame} returned by
#'   \code{\link{alpha_dominance}}, with row names corresponding to
#'   sample IDs and a column named \code{"Gini"}. If a column named
#'   \code{"scarcity"} is present, it will be included in the output.
#'
#' @return A \code{data.frame} with one row per sample and columns:
#'   \itemize{
#'     \item \code{sample} — sample ID
#'     \item \code{Axis} — 1D ordination coordinate in [0, 1]
#'     \item \code{Gini} — dominance Gini for the sample
#'     \item \code{scarcity} — (optional) scarcity from \code{alpha_dom},
#'           if available
#'   }
#'
#' @details
#' Let \eqn{s_{\min}} and \eqn{s_{\max}} be the samples with minimum and
#' maximum Gini dominance, and let \eqn{d(s_{\min}, s_{\max})} be their
#' pairwise distance. For any sample \eqn{s}, its coordinate is defined as
#' \deqn{
#'   \text{Axis}(s) = \frac{d(s_{\min}, s)}{d(s_{\min}, s_{\max})},
#' }
#' with \eqn{\text{Axis}(s_{\min}) = 0} and \eqn{\text{Axis}(s_{\max}) = 1}.
#'
#' If all valid Gini values are identical (i.e., no minimum/maximum
#' separation), or if the distance between the two poles is zero, the
#' function stops with an error.
#'
#' @examples
#' set.seed(1)
#' # fake distance matrix
#' D <- as.matrix(dist(matrix(rnorm(50), nrow = 10)))
#' rownames(D) <- colnames(D) <- paste0("Sample_", 1:10)
#'
#' # fake alpha_dominance result
#' alpha_res <- data.frame(
#'   scarcity = runif(10),
#'   Gini     = runif(10),
#'   row.names = rownames(D)
#' )
#'
#' polar_ordination(D, alpha_res)
#'
#' @export
polar_ordination <- function(dist_obj, alpha_dom) {
  # Convert dist to matrix if needed
  if (inherits(dist_obj, "dist")) {
    dist_mat <- as.matrix(dist_obj)
  } else {
    dist_mat <- as.matrix(dist_obj)
  }
  
  # Basic checks
  if (is.null(rownames(dist_mat)) || is.null(colnames(dist_mat))) {
    stop("dist_obj must have row and column names (sample IDs).")
  }
  alpha_dom <- as.data.frame(alpha_dom)
  if (is.null(rownames(alpha_dom))) {
    stop("alpha_dom must have row names corresponding to sample IDs.")
  }
  if (!("Gini" %in% colnames(alpha_dom))) {
    stop("alpha_dom must contain a column named 'Gini'.")
  }
  
  # Align samples
  common_samples <- intersect(rownames(dist_mat), rownames(alpha_dom))
  if (length(common_samples) < 2) {
    stop("Fewer than 2 common samples between dist_obj and alpha_dom.")
  }
  
  dist_mat_sub <- dist_mat[common_samples, common_samples, drop = FALSE]
  alpha_sub    <- alpha_dom[common_samples, , drop = FALSE]
  
  # Extract Gini and handle NAs
  gini_vec <- alpha_sub$Gini
  names(gini_vec) <- common_samples
  valid <- !is.na(gini_vec)
  
  if (sum(valid) < 2) {
    stop("Fewer than 2 samples with non-NA Gini values.")
  }
  
  gini_vec <- gini_vec[valid]
  dist_mat_sub <- dist_mat_sub[names(gini_vec), names(gini_vec), drop = FALSE]
  alpha_sub    <- alpha_sub[names(gini_vec), , drop = FALSE]
  
  # Identify poles: min and max Gini
  pole_min <- names(gini_vec)[which.min(gini_vec)]
  pole_max <- names(gini_vec)[which.max(gini_vec)]
  
  if (identical(pole_min, pole_max)) {
    stop("All Gini values are identical; cannot define distinct poles.")
  }
  
  d_poles <- dist_mat_sub[pole_min, pole_max]
  if (d_poles == 0) {
    stop("The two poles have zero distance. Cannot define axis.")
  }
  
  # Ordination: distance to min-Gini pole, normalized by distance between poles
  ordination_values <- dist_mat_sub[pole_min, ] / d_poles
  ordination_values[pole_min] <- 0
  ordination_values[pole_max] <- 1
  
  ordination_values <- as.numeric(ordination_values)
  names(ordination_values) <- colnames(dist_mat_sub)
  
  # Build result
  res <- data.frame(
    sample = names(ordination_values),
    Axis   = ordination_values,
    Gini   = alpha_sub[names(ordination_values), "Gini"]
  )
  
  # Attach scarcity if present
  if ("scarcity" %in% colnames(alpha_sub)) {
    res$scarcity <- alpha_sub[names(ordination_values), "scarcity"]
  }
  
  rownames(res) <- res$sample
  res
}
