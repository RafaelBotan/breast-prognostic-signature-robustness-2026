# Data dictionary

All CSV files are aggregate or derived. None has one row per patient.

## Final revision outputs (`results/`)

| File | Rows | Unit and purpose |
|---|---:|---|
| `cohort_censoring.csv` | 6 | One cohort per row: endpoint, analyzed sample, events, censoring, and reverse Kaplan-Meier follow-up. |
| `critical_gene_summary.csv` | 171 | One signature-gene summary per row across evaluated cohorts. |
| `e2_transport_fixed_orientation.csv` | 56 | One signature-target-platform-normalization mode per row after fixing prognostic direction a priori. |
| `envelope_neff_24cells.csv` | 24 | One architecture-by-cohort cell with panel size, effective gene number, intact C-index, and gene-loss envelope. |
| `illustrative_cost_sensitivity.csv` | 18 | Hypothetical analytical-well scenarios; this is a sensitivity analysis, not an economic evaluation. |
| `leave_one_gene_out.csv` | 962 | One cohort-signature-gene omission per row with the resulting C-index decrease. |
| `matched_size_mc_draws.csv` | 27,000 | One attempted Monte Carlo gene subset per row; these are simulation draws, not patients or biological replicates. |
| `matched_size_mc_stability.csv` | 108 | Prefix summaries after 30, 100, 500, and 1,000 attempted draws. |
| `matched_size_mc_summary.csv` | 27 | One cohort-architecture-target-size summary with valid/nonestimable counts and empirical percentiles. |
| `neff_envelope_correlations.csv` | 7 | Overall and cohort-specific descriptive Spearman correlations without independence-based P values. |
| `orientation_audit.csv` | 24 | One architecture-by-cohort cell recording the prespecified prognostic direction used before perturbation. |

`results/e2_transport_fixed_orientation.csv` is byte-identical to `data/reference_aggregate/e2_transport.csv`. Fixing prognostic direction a priori did not change any of the 56 E2 cells; the revised analysis makes the orientation rule explicit while preserving the aggregate values.

## Aggregate reference inputs (`data/reference_aggregate/`)

| File | Rows | Purpose |
|---|---:|---|
| `e2_transport.csv` | 56 | Original aggregate transportability output used as a deterministic reference. |
| `stress_curves.csv` | 2,255 | Aggregate degradation-curve points by cohort, architecture, and gene-loss step. |
| `v2_envelope_abs_neff.csv` | 20 | Original exploratory envelope and effective-gene-number summaries. |
| `v2_matched_size.csv` | 10 | Original aggregate 24-gene matched-size summaries. |
| `v2b_matched_size15.csv` | 15 | Original aggregate 15-gene matched-size summaries. |
| `v4_heldout_GSE7390.csv` | 4 | Held-out cohort aggregate results. |
| `v4_heldout_matched.csv` | 2 | Held-out cohort aggregate matched-size summaries. |

## Key variables

- `c_full`: Harrell concordance index for the intact harmonized score.
- `envelope`: number of genes removable in decreasing order of absolute published weight before Harrell C decreases by at least 0.02; at least three genes must remain.
- `n_eff`: effective gene number, `(sum(abs(w))^2) / sum(w^2)`, which decreases as weight becomes concentrated in fewer genes.
- `c_drop`: intact Harrell C minus Harrell C after omitting one gene.
- `direction_reversed`: whether the published score direction required a sign reversal fixed before perturbation analysis.
- `nonestimable`: attempted Monte Carlo subset that did not satisfy the prespecified concordance estimability checks.
- `frozen_within`: external means and standard deviations estimated from other cohorts in the same broad technology family, never from the target cohort.
- `frozen_cross`: fixed SCAN-B reference; for microarray targets this is a cross-technology stress, not an isolated platform effect.
