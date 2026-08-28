# Robustness of Harmonized Multigene Prognostic Scores to Gene Loss and Single-Sample Scoring Across Six Breast Cancer Cohorts

Reproducibility package for the major revision of manuscript **CCI-26-00256**, submitted to *JCO Clinical Cancer Informatics*.

## Scope

This repository contains the code and aggregate or derived outputs needed to audit the revised analyses. It evaluates four harmonized, linear breast-cancer prognostic score architectures under structured gene loss and external-reference normalization. The implementations are research approximations; they are not the proprietary commercial assays and must not be used to modify a validated clinical test.

No patient-level expression, survival, identifier, or report data are included. The underlying cohorts are public, but their processed patient-level matrices must be obtained from the original repositories and prepared locally.

## Authors

- Rafael de Negreiros Botan — Postgraduate Program in Medical Sciences, School of Medicine, Universidade de Brasilia; ORCID 0000-0002-7290-5824
- Joao Batista de Sousa — Faculty of Medicine, Universidade de Brasilia

## Repository map

```text
code/original/            Original analysis scripts retained as executable provenance
code/revision/            Revised and additional analyses requested during peer review
data/reference_aggregate/ Aggregate outputs used as deterministic reference inputs
results/                  Final aggregate or derived outputs reported in the revision
figures/                  Figures 1 and 2 in PNG, TIFF, SVG, and PDF
environment/              R version and installed-package session information
docs/                     Reproducibility quality checks
DATA_DICTIONARY.md        File-level data dictionary and key variable definitions
MANIFEST.sha256           SHA-256 custody manifest for every published file
```

Internal editorial strategy memoranda, reviewer correspondence, the manuscript, and patient-level observations are deliberately excluded. Three submitted-versus-revised comparison tables used during internal quality control are also excluded; the repository publishes the final analysis state rather than editorial working files.

## Public cohort sources

| Cohort | Source |
|---|---|
| SCAN-B | GEO GSE96058 |
| TCGA-BRCA | GDC or cBioPortal |
| METABRIC | cBioPortal or EGA EGAS00000000083 |
| GSE20685 | GEO GSE20685 |
| GSE1456 | GEO GSE1456 |
| GSE7390 | GEO GSE7390 |

The revision scripts expect one local folder per cohort (`SCANB`, `TCGA_BRCA`, `METABRIC`, `GSE20685`, `GSE1456`, and `GSE7390`) containing the processed public-data files `expression_genelevel_preZ.parquet` and `clinical_FINAL.parquet`. These files are not redistributed here.

## Reproduce the revision

The final build used **R 4.5.3**. Exact installed-package versions are recorded in `environment/sessionInfo.txt`. Principal packages include `arrow`, `survival`, `glmnet`, `genefu`, `ggplot2`, `patchwork`, `scales`, and `svglite`.

On Windows PowerShell:

```powershell
.\reproduce.ps1 -CohortDir 'D:\path\to\processed_public_cohorts'
```

On Bash:

```bash
COHORT_DIR=/path/to/processed_public_cohorts ./reproduce.sh
```

Outputs are written under `outputs/`. The Monte Carlo analysis uses seed `20260827` and 1,000 attempted gene subsets per architecture, cohort, and target size. Prognostic direction is fixed before perturbation or outcome evaluation. The summaries report empirical 10th and 90th percentiles and do not treat Monte Carlo draws as independent biological replicates.

The published `results/` directory contains the verified reference outputs. `docs/REPRODUCIBILITY_QA.md` records the deterministic checks performed before publication.

## Data and code availability

Only code and aggregate or derived data are released. No individual patient data are shared. See `CODE_DATA_AVAILABILITY.md` for manuscript-ready wording.

## License and citation

- Code: MIT (`LICENSE`).
- Aggregate data and documentation: CC BY 4.0 (`LICENSE-DATA.md`).
- Citation metadata: `CITATION.cff`.

