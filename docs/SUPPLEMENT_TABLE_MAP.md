# Supplement-to-repository evidence map

This map identifies the aggregate evidence and executable code supporting each table and figure in the revised Data Supplement. All listed CSV rows are cohort-, architecture-, gene-, normalization-, or simulation-level; none is one row per patient.

| Supplement item | Primary aggregate evidence | Principal code |
|---|---|---|
| Table S1: cohorts, end points, events, censoring, follow-up | `results/cohort_censoring.csv` | `code/revision/reviewer1_revision_analysis.R` |
| Table S2: intact Harrell/Uno concordance | `data/reference_aggregate/v3_uno_vs_harrell.csv`; `data/reference_aggregate/v4_heldout_GSE7390.csv` | `code/original/stress_v3_uno.R`; `code/original/stress_v4_heldout_and_fig.R` |
| Table S3: genes mapped and present | `data/reference_aggregate/v2_envelope_abs_neff.csv`; `data/reference_aggregate/v4_heldout_GSE7390.csv` | `code/original/stress_v2_matched_size.R`; `code/original/stress_v4_heldout_and_fig.R` |
| Table S4: primary envelope and bootstrap interval | `data/reference_aggregate/v2_bootstrap_env.csv`; `data/reference_aggregate/v4_heldout_GSE7390.csv` | `code/original/stress_v2_matched_size.R`; `code/original/stress_v4_heldout_and_fig.R` |
| Table S5: three removal sequences | `data/reference_aggregate/stress_envelope.csv` | `code/original/locked_stress_test.R` |
| Table S6: thresholds 0.01, 0.02, and 0.03 | `data/reference_aggregate/stress_curves.csv` | `code/original/locked_stress_test.R` |
| Table S7: matched-size analysis | `results/matched_size_mc_summary.csv` | `code/revision/reviewer1_revision_analysis.R`; `code/revision/postprocess_matched_outputs.R` |
| Table S8: matched-size stability and revision reconciliation | `results/matched_size_mc_stability.csv`; `results/matched_revision_comparison.csv` | `code/revision/postprocess_matched_outputs.R` |
| Table S9: Harrell/Uno envelope comparison | `data/reference_aggregate/v3_uno_vs_harrell.csv`; `data/reference_aggregate/v4_heldout_GSE7390.csv` | `code/original/stress_v3_uno.R`; `code/original/stress_v4_heldout_and_fig.R` |
| Table S10: effective gene number and envelope | `results/neff_envelope_correlations.csv`; `results/envelope_neff_24cells.csv` | `code/revision/reviewer1_revision_analysis.R` |
| Table S11: leave-one-gene-out results | `results/leave_one_gene_out.csv`; `results/critical_gene_summary.csv` | `code/revision/reviewer1_revision_analysis.R` |
| Table S12: frozen-reference scoring | `results/e2_transport_fixed_orientation.csv`; `results/e2_submitted_reconciliation.csv` | `code/revision/stress_e2_platform_revision.R` |
| Table S13: analytical-well sensitivity | `results/illustrative_cost_sensitivity.csv` | `code/revision/reviewer1_revision_analysis.R` |
| Table S14: exploratory rederivation | `data/reference_aggregate/rederive_results.csv` | `code/original/rederive_transportable.R` |
| Table S15: submitted-versus-regenerated reconciliation | `results/submitted_result_reconciliation.csv`; `results/e2_submitted_reconciliation.csv`; `results/matched_revision_comparison.csv` | `code/revision/reviewer1_revision_analysis.R`; `code/revision/stress_e2_platform_revision.R`; `code/revision/postprocess_matched_outputs.R` |
| Figure 1 | `data/reference_aggregate/stress_curves.csv` | `code/revision/make_figure1_revision.R` |
| Figure 2 | `results/matched_size_mc_summary.csv`; `results/envelope_neff_24cells.csv` | `code/revision/make_figure2_revision.R` |

The scripts under `code/original/` preserve submitted-analysis provenance. The current revised matched-size and fixed-orientation analyses are under `code/revision/`; when values changed during peer review, the revised values in `results/` are authoritative.
