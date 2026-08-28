suppressMessages({library(ggplot2); library(patchwork)})

project_dir <- Sys.getenv(
  "STRESS_PROJECT_DIR",
  normalizePath(".", winslash = "/", mustWork = TRUE)
)
data_dir <- Sys.getenv(
  "REVISION_OUTPUT_DIR",
  file.path(project_dir, "revision_2026-08-27", "analysis", "outputs")
)
out <- Sys.getenv(
  "FIGURE_OUTPUT_DIR",
  file.path(project_dir, "revision_2026-08-27", "figures")
)
dir.create(out, recursive=TRUE, showWarnings=FALSE)

matched <- read.csv(file.path(data_dir, "matched_size_mc_summary.csv"), stringsAsFactors=FALSE)
cells <- read.csv(file.path(data_dir, "envelope_neff_24cells.csv"), stringsAsFactors=FALSE)

architecture_labels <- c(
  GGI="GGI (uniform)",
  MammaPrint="70-gene-like (near-uniform)",
  CorePAM="24-gene (concentrated)",
  OncotypeDX="21-gene-like (concentrated)"
)
palette <- c(
  "GGI (uniform)"="#0072B2",
  "70-gene-like (near-uniform)"="#56B4E9",
  "24-gene (concentrated)"="#D55E00",
  "21-gene-like (concentrated)"="#A33B1E"
)
cohort_order <- c("SCANB", "TCGA_BRCA", "METABRIC", "GSE20685", "GSE1456", "GSE7390")

large24 <- matched[matched$target_genes==24 & matched$signature %in% c("GGI", "MammaPrint"), ]
native24 <- cells[cells$signature=="CorePAM", c("cohort", "signature", "envelope")]
names(native24)[names(native24)=="envelope"] <- "median_envelope"
native24$p10 <- native24$median_envelope
native24$p90 <- native24$median_envelope
plot_a <- rbind(
  large24[, c("cohort", "signature", "median_envelope", "p10", "p90")],
  native24[, c("cohort", "signature", "median_envelope", "p10", "p90")]
)
plot_a$Architecture <- factor(architecture_labels[plot_a$signature], levels=architecture_labels[c("CorePAM", "MammaPrint", "GGI")])
plot_a$cohort_index <- match(plot_a$cohort, cohort_order)
plot_a$y_base <- c(
  "24-gene (concentrated)"=1,
  "70-gene-like (near-uniform)"=2,
  "GGI (uniform)"=3
)[as.character(plot_a$Architecture)]
plot_a$y <- plot_a$y_base + seq(-0.24, 0.24, length.out=6)[plot_a$cohort_index]
plot_a$Status <- ifelse(plot_a$cohort=="GSE7390", "Held-out cohort", "Exploratory cohort")

pA <- ggplot(plot_a, aes(color=Architecture)) +
  geom_segment(
    data=plot_a[plot_a$signature!="CorePAM", ],
    aes(x=p10, xend=p90, y=y, yend=y), linewidth=0.7, alpha=0.75
  ) +
  geom_point(aes(x=median_envelope, y=y, shape=Status), size=2.1, stroke=0.8) +
  scale_color_manual(values=palette) +
  scale_shape_manual(values=c("Exploratory cohort"=16, "Held-out cohort"=1)) +
  scale_x_continuous(breaks=seq(0, 21, 3), limits=c(-0.5, 21.5), expand=c(0,0)) +
  scale_y_continuous(
    breaks=1:3,
    labels=c("24-gene\n(concentrated, native)", "70-gene-like\n(reduced to 24)", "GGI\n(reduced to 24)"),
    limits=c(0.55, 3.45)
  ) +
  labs(
    title="A  Matched panel size (24 genes)",
    subtitle="Median (P10-P90), 1,000 draws per cohort",
    x="Genes removable before C-index decreases by 0.02", y=NULL, shape=NULL
  ) +
  theme_classic(base_size=10, base_family="sans") +
  theme(
    legend.position="bottom",
    legend.title=element_blank(),
    axis.text.y=element_text(size=8),
    plot.subtitle=element_text(size=8, color="#444444"),
    plot.margin=margin(8, 12, 6, 8)
  ) +
  guides(color="none", shape=guide_legend(override.aes=list(color="#333333")))

cells$Architecture <- factor(architecture_labels[cells$signature], levels=architecture_labels)
pB <- ggplot(cells, aes(n_eff, envelope, color=Architecture, shape=Architecture)) +
  geom_point(size=2.4, alpha=0.90) +
  annotate("text", x=10, y=84, hjust=0, label="Spearman rho = 0.831", size=3.1, color="#222222") +
  scale_color_manual(values=palette) +
  scale_x_continuous(breaks=c(10, 20, 40, 60, 80), limits=c(7, 92)) +
  scale_y_continuous(breaks=seq(0, 90, 15), limits=c(-2, 90)) +
  labs(
    title="B  Effective gene number and envelope",
    subtitle="Each point is one architecture in one cohort (24 dependent cells)",
    x="Effective gene number (n_eff)",
    y="Gene-loss envelope (genes removable)",
    color=NULL, shape=NULL
  ) +
  theme_classic(base_size=10, base_family="sans") +
  theme(
    legend.position="right",
    legend.text=element_text(size=7.5),
    legend.key.height=grid::unit(1.0, "lines"),
    plot.subtitle=element_text(size=8, color="#444444"),
    plot.margin=margin(8, 8, 6, 12)
  )

p <- pA + pB + plot_layout(widths=c(1.08, 1.42))
ggsave(file.path(out, "Figure2_revision.pdf"), p, width=10.8, height=4.5, device=cairo_pdf)
ggsave(file.path(out, "Figure2_revision.svg"), p, width=10.8, height=4.5, device=svg)
ggsave(file.path(out, "Figure2_revision.tiff"), p, width=10.8, height=4.5, dpi=400, compression="lzw")
ggsave(file.path(out, "Figure2_revision.png"), p, width=10.8, height=4.5, dpi=300, bg="white")
cat("Figure 2 complete; panel A rows:", nrow(plot_a), "; panel B rows:", nrow(cells), "\n")
