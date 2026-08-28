# Figure 1 (publication quality, English, hi-res) — gene-loss degradation curves.
suppressMessages(library(ggplot2))
project_dir <- Sys.getenv(
  "STRESS_PROJECT_DIR",
  normalizePath(".", winslash = "/", mustWork = TRUE)
)
aggregate_dir <- Sys.getenv(
  "SUBMITTED_AGGREGATE_DIR",
  dirname(project_dir)
)
OUT <- Sys.getenv("FIGURE_OUTPUT_DIR", project_dir)
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
cur <- read.csv(file.path(aggregate_dir, "stress_curves.csv"), check.names=FALSE)

d <- cur[cur$loss %in% c("weight","full"), ]
# rotulos do manuscrito
relab <- c(CorePAM="24-gene (concentrated)", OncotypeDX="21-gene-like (concentrated)",
           GGI="GGI (uniform)", MammaPrint="70-gene-like (near-uniform)")
d$Signature <- factor(relab[d$sig], levels=c("GGI (uniform)","70-gene-like (near-uniform)","24-gene (concentrated)","21-gene-like (concentrated)"))
cohlab <- c(SCANB="SCAN-B", TCGA_BRCA="TCGA-BRCA", METABRIC="METABRIC", GSE20685="GSE20685", GSE1456="GSE1456")
d$Cohort <- factor(cohlab[d$cohort], levels=c("SCAN-B","TCGA-BRCA","METABRIC","GSE20685","GSE1456"))

pal <- c("GGI (uniform)"="#1b6ca8","70-gene-like (near-uniform)"="#4ea1d3",
         "24-gene (concentrated)"="#d95f02","21-gene-like (concentrated)"="#a33b1e")

p <- ggplot(d, aes(frac_retained, cindex, color=Signature)) +
  geom_line(linewidth=0.7) + geom_point(size=0.6) +
  facet_wrap(~Cohort, scales="free_y", nrow=2) +
  scale_x_reverse(labels=scales::percent) +
  scale_color_manual(values=pal) +
  labs(x="Fraction of genes retained (highest-weight genes removed first)",
       y="Harrell C-index", color="Architecture") +
  theme_bw(base_size=11) +
  theme(legend.position="bottom", panel.grid.minor=element_blank(),
        strip.background=element_rect(fill="grey92", color=NA),
        legend.key.width=unit(1.1,"lines"))

ggsave(file.path(OUT,"Figure1.pdf"), p, width=9, height=6, device=cairo_pdf)
ggsave(file.path(OUT,"Figure1.tiff"), p, width=9, height=6, dpi=400, compression="lzw")
ggsave(file.path(OUT,"Figure1.png"),  p, width=9, height=6, dpi=300, bg="white")
cat("OK Figure1 pdf/tiff/png salvos\n")
