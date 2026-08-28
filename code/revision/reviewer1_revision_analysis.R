# Additional analyses for the major revision of CCI-26-00256.
#
# Purpose:
#   1. describe censoring by cohort;
#   2. replace the 30-draw language with a reproducible Monte Carlo analysis
#      (1,000 gene-set draws; descriptive, not biological replicates);
#   3. quantify the association between effective gene number and envelope;
#   4. identify genes whose single omission changes Harrell's C by >= 0.02;
#   5. provide a transparent, illustrative cost/QC sensitivity calculation.
#
# Public, de-identified cohorts only. Patient-level data remain local and are
# never written. Every output is aggregated by cohort, signature, gene, or
# hypothetical assay scenario.

suppressPackageStartupMessages({
  library(arrow)
  library(GEOquery)
  library(genefu)
  library(survival)
})

set.seed(20260827)

cohort_dir <- Sys.getenv("COHORT_DIR", "")
if (!nzchar(cohort_dir)) {
  stop("Set COHORT_DIR to the processed public-cohort directory.")
}
project_dir <- Sys.getenv(
  "STRESS_PROJECT_DIR",
  normalizePath(".", winslash = "/", mustWork = TRUE)
)
analysis_dir <- file.path(project_dir, "revision_2026-08-27", "analysis")
out_dir <- Sys.getenv(
  "REVISION_OUTPUT_DIR",
  file.path(analysis_dir, "outputs")
)
cache_dir <- Sys.getenv(
  "REVISION_CACHE_DIR",
  file.path(analysis_dir, "cache_geo")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

B_MATCH <- as.integer(Sys.getenv("B_MATCH", "1000"))
THRESHOLD <- 0.02

log_line <- function(...) {
  cat(sprintf(...), "\n")
  flush.console()
}

# ---- Signature definitions: identical to the submitted analysis ----
corepam <- c(
  EXO1 = 0.2174, NAT1 = -0.2052, BLVRA = -0.1848,
  ACTR3B = -0.1405, MIA = -0.1193, MYBL2 = -0.1183,
  PTTG1 = 0.1140, MDM2 = -0.1074, SFRP1 = -0.0907,
  GPR160 = -0.0841, FOXC1 = 0.0717, PHGDH = 0.0620,
  MYC = -0.0424, ESR1 = 0.0396, KRT5 = -0.0277,
  CXXC5 = 0.0266, PGR = -0.0262, KRT17 = -0.0163,
  CENPF = 0.0150, FGFR4 = 0.0110, BCL2 = -0.0105,
  MLPH = -0.0099, ERBB2 = -0.0082, GRB7 = -0.0053
)

data(sig.oncotypedx)
data(sig.ggi)
data(sig.gene70)

od <- sig.oncotypedx
od <- od[od$group != "reference", ]
oncotype <- setNames(
  od$weight * ifelse(od$group == "estrogen", -1, 1),
  od$symbol
)

gg <- sig.ggi[!is.na(sig.ggi$HUGO.gene.symbol), ]
ggi <- setNames(
  ifelse(gg$grade == max(gg$grade, na.rm = TRUE), 1, -1),
  gg$HUGO.gene.symbol
)
ggi <- ggi[!duplicated(names(ggi))]

g70 <- sig.gene70[!is.na(sig.gene70$HUGO.gene.symbol), ]
mammaprint <- setNames(-g70$correlation, g70$HUGO.gene.symbol)
mammaprint <- mammaprint[!duplicated(names(mammaprint))]

sigs <- list(
  GGI = ggi,
  MammaPrint = mammaprint,
  CorePAM = corepam,
  OncotypeDX = oncotype
)

signature_labels <- c(
  GGI = "GGI (uniform)",
  MammaPrint = "70-gene-like (near-uniform)",
  CorePAM = "24-gene (concentrated)",
  OncotypeDX = "21-gene-like (concentrated)"
)

# survival::concordance treats larger predictors as longer survival. The score
# directions are fixed here from the published definitions, before any cohort
# is examined. The multiplier is therefore never chosen from the evaluation
# outcomes. This coding reproduces the submitted intact-score orientation in
# all 24 architecture-by-cohort cells.
concordance_multiplier <- c(
  GGI = -1,
  MammaPrint = 1,
  CorePAM = -1,
  OncotypeDX = -1
)

endpoint_labels <- c(
  SCANB = "overall survival",
  TCGA_BRCA = "overall survival",
  METABRIC = "overall survival",
  GSE20685 = "overall survival",
  GSE1456 = "overall survival",
  GSE7390 = "distant metastasis-free survival"
)

# ---- Shared functions ----
score_sig <- function(zmat, w) {
  present <- names(w)[names(w) %in% rownames(zmat)]
  if (length(present) < 3) return(rep(NA_real_, ncol(zmat)))
  wp <- w[present]
  as.numeric(
    colSums(zmat[present, , drop = FALSE] * wp) / sum(abs(wp))
  )
}

c_index <- function(time, event, score, timewt = "n") {
  ok <- is.finite(time) & is.finite(event) & is.finite(score)
  if (sum(ok) < 20 || length(unique(score[ok])) < 3) return(NA_real_)
  tryCatch(
    concordance(
      Surv(time[ok], event[ok]) ~ score[ok],
      timewt = timewt
    )$concordance,
    error = function(e) NA_real_
  )
}

effective_gene_number <- function(w) {
  a <- abs(w)
  if (sum(a) == 0) return(NA_real_)
  (sum(a)^2) / sum(a^2)
}

envelope_weight <- function(zmat, time, event, w, threshold = THRESHOLD) {
  # `w` must already have its direction fixed from the intact signature.
  # Never re-orient a sampled subset on the outcomes used to evaluate it.
  c_full <- c_index(time, event, score_sig(zmat, w))
  if (!is.finite(c_full)) return(c(c_full = NA_real_, envelope = NA_real_))
  n_genes <- length(w)
  removal_order <- names(sort(abs(w), decreasing = TRUE))
  max_removable <- n_genes - 3
  for (n_removed in seq_len(n_genes - 3)) {
    keep <- setdiff(names(w), removal_order[seq_len(n_removed)])
    c_reduced <- c_index(time, event, score_sig(zmat, w[keep]))
    if (is.finite(c_reduced) && (c_full - c_reduced) >= threshold) {
      max_removable <- n_removed - 1
      break
    }
  }
  c(c_full = c_full, envelope = max_removable)
}

reverse_km_median <- function(time, event) {
  ok <- is.finite(time) & is.finite(event) & time >= 0
  if (sum(ok) < 2) return(NA_real_)
  fit <- tryCatch(
    survfit(Surv(time[ok], 1 - event[ok]) ~ 1),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NA_real_)
  tab <- summary(fit)$table
  if (is.null(tab) || !("median" %in% names(tab))) return(NA_real_)
  unname(tab[["median"]])
}

load_local_cohort <- function(cohort) {
  expression_path <- file.path(
    cohort_dir, cohort, "expression_genelevel_preZ.parquet"
  )
  clinical_path <- file.path(cohort_dir, cohort, "clinical_FINAL.parquet")
  stopifnot(file.exists(expression_path), file.exists(clinical_path))

  ex <- as.data.frame(read_parquet(expression_path))
  genes <- ex$gene
  ex$gene <- NULL
  ex <- as.matrix(ex)
  rownames(ex) <- genes
  ex <- ex[, !grepl("REPL$", colnames(ex)), drop = FALSE]

  clinical <- as.data.frame(read_parquet(clinical_path))
  key <- if (
    sum(colnames(ex) %in% clinical$patient_id) >
      sum(colnames(ex) %in% clinical$sample_id)
  ) "patient_id" else "sample_id"
  matched <- match(colnames(ex), clinical[[key]])
  time <- clinical$os_time_months[matched]
  event <- clinical$os_event[matched]

  zmat <- t(scale(t(ex)))
  zmat <- zmat[is.finite(rowSums(zmat)), , drop = FALSE]
  list(
    cohort = cohort,
    z = zmat,
    time = time,
    event = event,
    n_expression = ncol(ex)
  )
}

load_gse7390 <- function() {
  log_line("Loading GSE7390 from GEO (cached under revision analysis)...")
  gse <- getGEO(
    "GSE7390",
    GSEMatrix = TRUE,
    AnnotGPL = TRUE,
    destdir = cache_dir
  )[[1]]
  ex <- exprs(gse)
  feature <- fData(gse)
  pheno <- pData(gse)

  symbol_col <- colnames(feature)[
    grep("symbol", colnames(feature), ignore.case = TRUE)
  ][1]
  if (is.na(symbol_col)) stop("No gene-symbol column found for GSE7390")

  find_col <- function(pattern) {
    hits <- colnames(pheno)[grep(pattern, colnames(pheno), ignore.case = TRUE)]
    if (length(hits)) hits[1] else NA_character_
  }
  parse_numeric <- function(x) {
    as.numeric(gsub(".*: *", "", as.character(x)))
  }

  time_col <- find_col("t.dmfs|t_dmfs|dmfs.time|time.dmfs")
  event_col <- find_col("e.dmfs|e_dmfs|dmfs.event|event.dmfs")
  if (is.na(time_col) || is.na(event_col)) {
    characteristic_cols <- grep(
      "characteristics", colnames(pheno), ignore.case = TRUE
    )
    characteristic <- pheno[, characteristic_cols, drop = FALSE]
    flat <- apply(characteristic, 2, as.character)
    time_hits <- which(apply(
      flat, 2,
      function(x) any(grepl("t.dmfs|tdm", x, ignore.case = TRUE))
    ))
    event_hits <- which(apply(
      flat, 2,
      function(x) any(grepl("e.dmfs|edm", x, ignore.case = TRUE))
    ))
    if (!length(time_hits) || !length(event_hits)) {
      stop("GSE7390 DMFS columns not found")
    }
    # GEO reports t.dmfs in days. Convert only for readable follow-up output;
    # concordance itself is invariant to this positive rescaling of time.
    time <- parse_numeric(characteristic[, time_hits[1]]) / 30.4375
    event <- parse_numeric(characteristic[, event_hits[1]])
  } else {
    time <- parse_numeric(pheno[[time_col]]) / 30.4375
    event <- parse_numeric(pheno[[event_col]])
  }

  symbols <- sub(" ///.*", "", feature[[symbol_col]])
  keep <- symbols != "" & !is.na(symbols)
  ex <- ex[keep, , drop = FALSE]
  symbols <- symbols[keep]
  order_by_variance <- order(
    apply(ex, 1, var, na.rm = TRUE), decreasing = TRUE
  )
  ex <- ex[order_by_variance, , drop = FALSE]
  symbols <- symbols[order_by_variance]
  first_symbol <- !duplicated(symbols)
  ex <- ex[first_symbol, , drop = FALSE]
  rownames(ex) <- symbols[first_symbol]
  zmat <- t(scale(t(ex)))
  zmat <- zmat[is.finite(rowSums(zmat)), , drop = FALSE]

  list(
    cohort = "GSE7390",
    z = zmat,
    time = time,
    event = event,
    n_expression = ncol(ex)
  )
}

# ---- Run cohort-level analyses ----
cohort_order <- c(
  "SCANB", "TCGA_BRCA", "METABRIC", "GSE20685", "GSE1456", "GSE7390"
)

censoring_rows <- list()
orientation_rows <- list()
envelope_rows <- list()
critical_gene_rows <- list()
matched_draw_rows <- list()

write_checkpoint <- function(rows, filename) {
  if (!length(rows)) return(invisible(NULL))
  write.csv(
    do.call(rbind, rows),
    file.path(out_dir, filename),
    row.names = FALSE
  )
}

for (cohort in cohort_order) {
  dat <- if (cohort == "GSE7390") {
    load_gse7390()
  } else {
    load_local_cohort(cohort)
  }

  valid <- is.finite(dat$time) & is.finite(dat$event) & dat$time >= 0
  n_analyzed <- sum(valid)
  n_events <- sum(dat$event[valid] == 1)
  n_censored <- n_analyzed - n_events
  censoring_rows[[length(censoring_rows) + 1]] <- data.frame(
    cohort = cohort,
    endpoint = endpoint_labels[[cohort]],
    n_expression = dat$n_expression,
    n_analyzed = n_analyzed,
    events = n_events,
    censored = n_censored,
    censored_percent = round(100 * n_censored / n_analyzed, 1),
    median_followup_months_reverse_km = round(
      reverse_km_median(dat$time, dat$event), 1
    ),
    stringsAsFactors = FALSE
  )

  log_line(
    "== %s: n=%d, events=%d, censored=%d (%.1f%%) ==",
    cohort, n_analyzed, n_events, n_censored,
    100 * n_censored / n_analyzed
  )

  for (signature in names(sigs)) {
    w0 <- sigs[[signature]]
    present <- names(w0)[names(w0) %in% rownames(dat$z)]
    if (length(present) < 5) next
    w <- w0[present]

    raw_c <- c_index(dat$time, dat$event, score_sig(dat$z, w))
    oriented_w <- w * concordance_multiplier[[signature]]
    oriented_c <- c_index(
      dat$time, dat$event, score_sig(dat$z, oriented_w)
    )
    env <- envelope_weight(dat$z, dat$time, dat$event, oriented_w)

    orientation_rows[[length(orientation_rows) + 1]] <- data.frame(
      cohort = cohort,
      signature = signature,
      label = signature_labels[[signature]],
      n_genes = length(w),
      c_before_orientation = round(raw_c, 4),
      direction_reversed = concordance_multiplier[[signature]] == -1,
      c_after_orientation = round(oriented_c, 4),
      stringsAsFactors = FALSE
    )

    envelope_rows[[length(envelope_rows) + 1]] <- data.frame(
      cohort = cohort,
      signature = signature,
      label = signature_labels[[signature]],
      n_genes = length(w),
      n_eff = round(effective_gene_number(w), 4),
      c_full = round(env[["c_full"]], 4),
      envelope = as.integer(env[["envelope"]]),
      stringsAsFactors = FALSE
    )

    # Leave-one-gene-out analysis for an actionable QC example.
    for (gene in names(oriented_w)) {
      c_without <- c_index(
        dat$time,
        dat$event,
        score_sig(dat$z, oriented_w[names(oriented_w) != gene])
      )
      drop <- oriented_c - c_without
      critical_gene_rows[[length(critical_gene_rows) + 1]] <- data.frame(
        cohort = cohort,
        signature = signature,
        gene = gene,
        absolute_weight = abs(oriented_w[[gene]]),
        c_full = round(oriented_c, 4),
        c_without_gene = round(c_without, 4),
        c_drop = round(drop, 4),
        crosses_002 = is.finite(drop) && drop >= THRESHOLD,
        stringsAsFactors = FALSE
      )
    }

    # Conditional Monte Carlo draws. These quantify sensitivity to which genes
    # are selected; they are not independent patient-level or cohort-level
    # replicates and are not used for hypothesis tests.
    targets <- integer(0)
    if (signature %in% c("GGI", "MammaPrint")) targets <- c(24L, 15L)
    if (signature == "CorePAM") targets <- 15L
    if (cohort == "GSE7390") targets <- intersect(targets, 24L)

    for (target in targets) {
      if (length(w) <= target + 4) next
      for (draw in seq_len(B_MATCH)) {
        selected <- sample(names(oriented_w), target, replace = FALSE)
        draw_env <- envelope_weight(
          dat$z, dat$time, dat$event, oriented_w[selected]
        )
        matched_draw_rows[[length(matched_draw_rows) + 1]] <- data.frame(
          cohort = cohort,
          signature = signature,
          target_genes = target,
          draw = draw,
          c_full = round(draw_env[["c_full"]], 4),
          envelope = as.integer(draw_env[["envelope"]]),
          stringsAsFactors = FALSE
        )
      }
      log_line(
        "  %s matched to %d genes: %d Monte Carlo draws complete",
        signature, target, B_MATCH
      )
    }
  }
  # Durable cohort-level checkpoints: an interruption never discards completed
  # cohorts, and partial files are visibly labelled as such.
  write_checkpoint(censoring_rows, "partial_cohort_censoring.csv")
  write_checkpoint(orientation_rows, "partial_orientation_audit.csv")
  write_checkpoint(envelope_rows, "partial_envelope_neff.csv")
  write_checkpoint(critical_gene_rows, "partial_leave_one_gene_out.csv")
  write_checkpoint(matched_draw_rows, "partial_matched_size_mc_draws.csv")
  rm(dat)
  invisible(gc())
}

censoring <- do.call(rbind, censoring_rows)
orientation <- do.call(rbind, orientation_rows)
envelopes <- do.call(rbind, envelope_rows)
critical_genes <- do.call(rbind, critical_gene_rows)
matched_draws <- do.call(rbind, matched_draw_rows)

# ---- Reconcile the regenerated intact-score results to the submitted tables ----
source_dir <- Sys.getenv(
  "SUBMITTED_AGGREGATE_DIR",
  dirname(project_dir)
)
submitted_exploratory <- read.csv(
  file.path(source_dir, "v2_envelope_abs_neff.csv"),
  stringsAsFactors = FALSE
)
submitted_heldout <- read.csv(
  file.path(source_dir, "v4_heldout_GSE7390.csv"),
  stringsAsFactors = FALSE
)
submitted <- rbind(
  data.frame(
    cohort = submitted_exploratory$cohort,
    signature = submitted_exploratory$sig,
    submitted_c_full = submitted_exploratory$c_full,
    submitted_envelope = submitted_exploratory$env_abs,
    stringsAsFactors = FALSE
  ),
  data.frame(
    cohort = submitted_heldout$cohort,
    signature = submitted_heldout$sig,
    submitted_c_full = submitted_heldout$c_harrell,
    submitted_envelope = submitted_heldout$env_abs,
    stringsAsFactors = FALSE
  )
)
reconciliation <- merge(
  submitted,
  envelopes[, c("cohort", "signature", "c_full", "envelope")],
  by = c("cohort", "signature"),
  all = TRUE
)
reconciliation$c_difference <- round(
  reconciliation$c_full - reconciliation$submitted_c_full,
  6
)
reconciliation$envelope_difference <-
  reconciliation$envelope - reconciliation$submitted_envelope
reconciliation$matches_submitted <- with(
  reconciliation,
  is.finite(c_difference) & abs(c_difference) <= 0.0001 &
    is.finite(envelope_difference) & envelope_difference == 0
)
write.csv(
  reconciliation,
  file.path(out_dir, "submitted_result_reconciliation.csv"),
  row.names = FALSE
)
if (nrow(reconciliation) != 24 || !all(reconciliation$matches_submitted)) {
  stop(
    "Regenerated intact-score results do not reconcile to all 24 submitted cells; ",
    "inspect submitted_result_reconciliation.csv before interpretation."
  )
}

# ---- Correlation requested by Reviewer 1 ----
overall_rho <- suppressWarnings(cor(
  envelopes$n_eff,
  envelopes$envelope,
  method = "spearman",
  use = "complete.obs"
))

correlation_rows <- list(data.frame(
  scope = "all 24 architecture-by-cohort cells",
  n_cells = nrow(envelopes),
  spearman_rho = unname(overall_rho),
  interpretation = paste(
    "Descriptive only: cells share cohorts and architectures.",
    "The matched-size experiment, not a P value, addresses panel size."
  ),
  stringsAsFactors = FALSE
))

for (cohort in cohort_order) {
  x <- envelopes[envelopes$cohort == cohort, ]
  rho <- suppressWarnings(cor(
    x$n_eff, x$envelope, method = "spearman", use = "complete.obs"
  ))
  correlation_rows[[length(correlation_rows) + 1]] <- data.frame(
    scope = paste0("within ", cohort),
    n_cells = nrow(x),
    spearman_rho = unname(rho),
    interpretation = "Four architectures; descriptive cohort-level sensitivity.",
    stringsAsFactors = FALSE
  )
}
correlations <- do.call(rbind, correlation_rows)

# ---- Summaries and convergence of the Monte Carlo analysis ----
summarize_draws <- function(x) {
  target <- attr(x, "target")
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(c(
      median = NA_real_, p10 = NA_real_, p90 = NA_real_,
      mean = NA_real_, sd = NA_real_, ceiling_fraction = NA_real_
    ))
  }
  c(
    median = median(x),
    # The envelope is a discrete gene count. Use empirical quantiles so the
    # reported limits are observed counts rather than interpolated fractions.
    p10 = unname(quantile(x, 0.10, type = 1)),
    p90 = unname(quantile(x, 0.90, type = 1)),
    mean = mean(x),
    sd = sd(x),
    ceiling_fraction = mean(x == (target - 3))
  )
}

matched_summary_rows <- list()
stability_rows <- list()
group_keys <- unique(matched_draws[, c("cohort", "signature", "target_genes")])
for (i in seq_len(nrow(group_keys))) {
  key <- group_keys[i, ]
  x <- matched_draws[
    matched_draws$cohort == key$cohort &
      matched_draws$signature == key$signature &
      matched_draws$target_genes == key$target_genes,
  ]
  env_values <- x$envelope
  attr(env_values, "target") <- key$target_genes
  stats <- summarize_draws(env_values)
  matched_summary_rows[[length(matched_summary_rows) + 1]] <- data.frame(
    cohort = key$cohort,
    signature = key$signature,
    target_genes = key$target_genes,
    n_draws = nrow(x),
    median_envelope = stats[["median"]],
    p10 = stats[["p10"]],
    p90 = stats[["p90"]],
    mean_envelope = round(stats[["mean"]], 3),
    sd_envelope = round(stats[["sd"]], 3),
    ceiling_fraction = round(stats[["ceiling_fraction"]], 3),
    stringsAsFactors = FALSE
  )

  for (n_prefix in c(30L, 100L, 500L, B_MATCH)) {
    n_use <- min(n_prefix, nrow(x))
    prefix <- x$envelope[seq_len(n_use)]
    prefix <- prefix[is.finite(prefix)]
    stability_rows[[length(stability_rows) + 1]] <- data.frame(
      cohort = key$cohort,
      signature = key$signature,
      target_genes = key$target_genes,
      draws_used = n_use,
      median_envelope = if (length(prefix)) median(prefix) else NA_real_,
      p10 = if (length(prefix)) unname(quantile(prefix, 0.10, type = 1)) else NA_real_,
      p90 = if (length(prefix)) unname(quantile(prefix, 0.90, type = 1)) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
}
matched_summary <- do.call(rbind, matched_summary_rows)
matched_stability <- unique(do.call(rbind, stability_rows))

# ---- Gene-level QC summary ----
# Build a compact summary directly so every output column is scalar.
critical_keys <- unique(critical_genes[, c("signature", "gene")])
critical_summary_rows <- vector("list", nrow(critical_keys))
for (i in seq_len(nrow(critical_keys))) {
  key <- critical_keys[i, ]
  x <- critical_genes[
    critical_genes$signature == key$signature &
      critical_genes$gene == key$gene,
  ]
  critical_summary_rows[[i]] <- data.frame(
    signature = key$signature,
    gene = key$gene,
    cohorts_evaluated = nrow(x),
    cohorts_crossing_002 = sum(x$crosses_002, na.rm = TRUE),
    mean_c_drop = round(mean(x$c_drop, na.rm = TRUE), 4),
    max_c_drop = round(max(x$c_drop, na.rm = TRUE), 4),
    stringsAsFactors = FALSE
  )
}
critical_summary <- do.call(rbind, critical_summary_rows)
critical_summary <- critical_summary[
  order(
    critical_summary$signature,
    -critical_summary$cohorts_crossing_002,
    -critical_summary$max_c_drop
  ),
]

# ---- Illustrative cost/QC sensitivity model ----
# This is not a health-economic analysis. It shows the explicit assumptions a
# laboratory would need to supply: number of genes protected (k), per-well
# failure probability (p), independence of duplicate failures, and the extra
# number of wells. Costs are normalized to one unit per gene measurement.
cost_grid <- expand.grid(
  panel_genes = c(15L, 24L),
  protected_genes = c(1L, 3L, 5L),
  per_well_failure = c(0.01, 0.05, 0.10),
  stringsAsFactors = FALSE
)
cost_grid <- cost_grid[
  cost_grid$protected_genes <= cost_grid$panel_genes,
]
cost_grid$prob_any_protected_gene_unmeasured_single <- with(
  cost_grid,
  1 - (1 - per_well_failure)^protected_genes
)
cost_grid$prob_any_protected_gene_unmeasured_duplicate <- with(
  cost_grid,
  1 - (1 - per_well_failure^2)^protected_genes
)
cost_grid$relative_failure_reduction <- with(
  cost_grid,
  1 - prob_any_protected_gene_unmeasured_duplicate /
    prob_any_protected_gene_unmeasured_single
)
cost_grid$additional_wells <- cost_grid$protected_genes
cost_grid$analytical_well_increase_percent <- with(
  cost_grid,
  100 * additional_wells / panel_genes
)
numeric_cols <- vapply(cost_grid, is.numeric, logical(1))
cost_grid[numeric_cols] <- lapply(cost_grid[numeric_cols], round, 6)

# ---- Write only aggregated outputs ----
write.csv(censoring, file.path(out_dir, "cohort_censoring.csv"), row.names = FALSE)
write.csv(orientation, file.path(out_dir, "orientation_audit.csv"), row.names = FALSE)
write.csv(envelopes, file.path(out_dir, "envelope_neff_24cells.csv"), row.names = FALSE)
write.csv(correlations, file.path(out_dir, "neff_envelope_correlations.csv"), row.names = FALSE)
write.csv(critical_genes, file.path(out_dir, "leave_one_gene_out.csv"), row.names = FALSE)
write.csv(critical_summary, file.path(out_dir, "critical_gene_summary.csv"), row.names = FALSE)
write.csv(matched_draws, file.path(out_dir, "matched_size_mc_draws.csv"), row.names = FALSE)
write.csv(matched_summary, file.path(out_dir, "matched_size_mc_summary.csv"), row.names = FALSE)
write.csv(matched_stability, file.path(out_dir, "matched_size_mc_stability.csv"), row.names = FALSE)
write.csv(cost_grid, file.path(out_dir, "illustrative_cost_sensitivity.csv"), row.names = FALSE)

log_line("\nOutputs written to: %s", out_dir)
log_line(
  "Spearman rho across %d cells: %.3f (descriptive only; no P value)",
  nrow(envelopes),
  unname(overall_rho)
)
log_line("DONE")
