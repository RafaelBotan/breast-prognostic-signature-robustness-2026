# Reproducibility package quality checks

Final local checks on 2026-08-27:

- All 12 R scripts parsed successfully under R 4.5.3.
- All 24 intact architecture-by-cohort results and 24 gene-loss envelopes were regenerated exactly from the processed public cohorts.
- The aggregate-only postprocessing script reproduced the matched-size summary and stability outputs byte for byte.
- Figure scripts regenerated PNG, TIFF, and SVG versions of Figures 1 and 2 byte for byte.
- Both SVG files parsed as valid XML.
- Both PDF figures contain vector text and graphics without embedded raster panels.
- The package contains only code, simulation-level or aggregate public-cohort results, and figures; no patient-level expression or survival rows are included.
- No machine-specific absolute path remains in the published scripts.

Full regeneration requires the separately obtained processed public cohorts described in `README.md`. The verified reference outputs are included under `results/` so readers can inspect the reported values without receiving patient-level data.

