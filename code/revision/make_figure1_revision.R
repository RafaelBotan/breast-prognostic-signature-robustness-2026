suppressMessages(library(ggplot2))

project_dir <- Sys.getenv(
  "STRESS_PROJECT_DIR",
  normalizePath(".", winslash = "/", mustWork = TRUE)
)
aggregate_dir <- Sys.getenv(
  "SUBMITTED_AGGREGATE_DIR",
  dirname(project_dir)
)
out <- Sys.getenv(
  "FIGURE_OUTPUT_DIR",
  file.path(project_dir, "revision_2026-08-27", "figures")
)
dir.create(out, recursive=TRUE, showWarnings=FALSE)
cur <- read.csv(file.path(aggregate_dir, "stress_curves.csv"), check.names=FALSE)

d <- cur[cur$loss %in% c("weight", "full"), ]
labels <- c(
  GGI="GGI (uniform)",
  MammaPrint="70-gene-like (near-uniform)",
  CorePAM="24-gene (concentrated)",
  OncotypeDX="21-gene-like (concentrated)"
)
d$Architecture <- factor(labels[d$sig], levels=labels)
cohort_labels <- c(
  SCANB="SCAN-B", TCGA_BRCA="TCGA-BRCA", METABRIC="METABRIC",
  GSE20685="GSE20685", GSE1456="GSE1456"
)
d$Cohort <- factor(cohort_labels[d$cohort], levels=cohort_labels)

palette <- c(
  "GGI (uniform)"="#0072B2",
  "70-gene-like (near-uniform)"="#56B4E9",
  "24-gene (concentrated)"="#D55E00",
  "21-gene-like (concentrated)"="#A33B1E"
)

p <- ggplot(d, aes(frac_retained, cindex, color=Architecture)) +
  geom_line(linewidth=0.75) +
  geom_point(size=0.65) +
  facet_wrap(~Cohort, scales="free_y", nrow=2) +
  scale_x_reverse(labels=scales::percent) +
  scale_color_manual(values=palette) +
  labs(
    x="Genes retained (highest absolute weight removed first)",
    y="Harrell C-index", color="Architecture"
  ) +
  theme_bw(base_size=11, base_family="sans") +
  theme(
    legend.position="bottom",
    panel.grid.minor=element_blank(),
    strip.background=element_rect(fill="#F2F2F2", color=NA),
    legend.key.width=grid::unit(1.0, "lines"),
    plot.margin=margin(8, 10, 6, 8)
  )

ggsave(file.path(out, "Figure1_revision.pdf"), p, width=9, height=6, device=cairo_pdf)
# The base Windows SVG device fails for this faceted plot under R 4.5.3;
# svglite writes the same editable vector content reproducibly.
ggsave(
  file.path(out, "Figure1_revision.svg"), p, width=9, height=6,
  device=svglite::svglite
)
ggsave(file.path(out, "Figure1_revision.tiff"), p, width=9, height=6, dpi=400, compression="lzw")
ggsave(file.path(out, "Figure1_revision.png"), p, width=9, height=6, dpi=300, bg="white")
cat("Figure 1 complete\n")
