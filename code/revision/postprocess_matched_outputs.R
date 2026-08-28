# Post-process the completed Monte Carlo draws without re-reading patient data.
# This script is deterministic and writes only aggregate summaries.

project_dir <- Sys.getenv(
  "STRESS_PROJECT_DIR",
  normalizePath(".", winslash = "/", mustWork = TRUE)
)
source_dir <- Sys.getenv(
  "SUBMITTED_AGGREGATE_DIR",
  dirname(project_dir)
)
out_dir <- Sys.getenv(
  "REVISION_OUTPUT_DIR",
  file.path(project_dir, "revision_2026-08-27", "analysis", "outputs")
)

draws <- read.csv(
  file.path(out_dir, "matched_size_mc_draws.csv"),
  stringsAsFactors = FALSE
)
orientation <- read.csv(
  file.path(out_dir, "orientation_audit.csv"),
  stringsAsFactors = FALSE
)

# Confirm that the completed run used the same direction implied by the fixed,
# pre-outcome coding now written into reviewer1_revision_analysis.R.
expected_reversal <- c(
  GGI = TRUE,
  MammaPrint = FALSE,
  CorePAM = TRUE,
  OncotypeDX = TRUE
)
observed <- as.logical(orientation$direction_reversed)
expected <- unname(expected_reversal[orientation$signature])
stopifnot(length(observed) == 24L, all(observed == expected))

keys <- unique(draws[, c("cohort", "signature", "target_genes")])
summary_rows <- list()
stability_rows <- list()

for (i in seq_len(nrow(keys))) {
  key <- keys[i, ]
  x <- draws[
    draws$cohort == key$cohort &
      draws$signature == key$signature &
      draws$target_genes == key$target_genes,
  ]
  values <- x$envelope[is.finite(x$envelope)]
  if (length(values) < 0.95 * nrow(x)) {
    stop("More than 5% of Monte Carlo draws are non-estimable in ",
         key$cohort, "/", key$signature, "/", key$target_genes)
  }
  ceiling <- key$target_genes - 3
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    cohort = key$cohort,
    signature = key$signature,
    target_genes = key$target_genes,
    n_draws_attempted = nrow(x),
    n_draws_valid = length(values),
    n_draws_nonestimable = nrow(x) - length(values),
    median_envelope = median(values),
    # The envelope is a discrete gene count; report observed empirical
    # quantiles rather than interpolated fractional gene counts.
    p10 = unname(quantile(values, 0.10, type = 1)),
    p90 = unname(quantile(values, 0.90, type = 1)),
    mean_envelope = round(mean(values), 3),
    sd_envelope = round(sd(values), 3),
    ceiling_fraction = round(mean(values == ceiling), 3),
    stringsAsFactors = FALSE
  )
  for (n_prefix in c(30L, 100L, 500L, 1000L)) {
    prefix_attempts <- x$envelope[seq_len(min(n_prefix, nrow(x)))]
    prefix <- prefix_attempts[is.finite(prefix_attempts)]
    stability_rows[[length(stability_rows) + 1]] <- data.frame(
      cohort = key$cohort,
      signature = key$signature,
      target_genes = key$target_genes,
      draws_attempted = length(prefix_attempts),
      draws_valid = length(prefix),
      median_envelope = median(prefix),
      p10 = unname(quantile(prefix, 0.10, type = 1)),
      p90 = unname(quantile(prefix, 0.90, type = 1)),
      stringsAsFactors = FALSE
    )
  }
}

summary_new <- do.call(rbind, summary_rows)
stability <- do.call(rbind, stability_rows)

old24 <- read.csv(
  file.path(source_dir, "v2_matched_size.csv"),
  stringsAsFactors = FALSE
)
old15 <- read.csv(
  file.path(source_dir, "v2b_matched_size15.csv"),
  stringsAsFactors = FALSE
)
old_heldout <- read.csv(
  file.path(source_dir, "v4_heldout_matched.csv"),
  stringsAsFactors = FALSE
)

old <- rbind(
  data.frame(
    cohort = old24$cohort,
    signature = old24$sig,
    target_genes = 24L,
    old_draws = old24$n_draws,
    old_median = old24$env_abs_matched_med,
    old_p10 = old24$env_abs_matched_lo,
    old_p90 = old24$env_abs_matched_hi
  ),
  data.frame(
    cohort = old15$cohort,
    signature = old15$sig,
    target_genes = 15L,
    old_draws = old15$n_draws,
    old_median = old15$env_abs_matched_med,
    old_p10 = old15$env_abs_matched_lo,
    old_p90 = old15$env_abs_matched_hi
  ),
  data.frame(
    cohort = old_heldout$cohort,
    signature = old_heldout$sig,
    target_genes = 24L,
    old_draws = old_heldout$n_draws,
    old_median = old_heldout$env_matched_med,
    old_p10 = old_heldout$env_matched_lo,
    old_p90 = old_heldout$env_matched_hi
  )
)

comparison <- merge(
  old,
  summary_new,
  by = c("cohort", "signature", "target_genes"),
  all = TRUE
)
comparison$median_change <- comparison$median_envelope - comparison$old_median
comparison$p10_change <- comparison$p10 - comparison$old_p10
comparison$p90_change <- comparison$p90 - comparison$old_p90

write.csv(
  summary_new,
  file.path(out_dir, "matched_size_mc_summary.csv"),
  row.names = FALSE
)
write.csv(
  stability,
  file.path(out_dir, "matched_size_mc_stability.csv"),
  row.names = FALSE
)
write.csv(
  comparison,
  file.path(out_dir, "matched_revision_comparison.csv"),
  row.names = FALSE
)

cat("postprocess complete:", nrow(summary_new), "cells\n")
cat(
  "median changed in",
  sum(comparison$median_change != 0, na.rm = TRUE),
  "of", nrow(comparison), "cells\n"
)
