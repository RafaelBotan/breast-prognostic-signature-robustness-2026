# Figura 2 do manuscrito v2 (JCO CCI) — corrigida: eixo y LINEAR para que os
# envelopes iguais a ZERO apareçam (a versao log10 os descartava silenciosamente).
suppressMessages({library(ggplot2); library(patchwork)})
# Portable path adaptation for the public repository; analytic logic is unchanged.
OUT <- Sys.getenv("REVISION_OUTPUT_DIR", "results")
A <- read.csv(file.path(OUT,"v2_envelope_abs_neff.csv"), stringsAsFactors=FALSE)
M <- read.csv(file.path(OUT,"v2_matched_size.csv"),      stringsAsFactors=FALSE)
H <- read.csv(file.path(OUT,"v4_heldout_GSE7390.csv"),   stringsAsFactors=FALSE)
HM<- read.csv(file.path(OUT,"v4_heldout_matched.csv"),   stringsAsFactors=FALSE)

A <- rbind(A[,c("cohort","sig","n_genes","env_abs","gini","n_eff")],
           H[,c("cohort","sig","n_genes","env_abs","gini","n_eff")])
M <- rbind(M[,c("cohort","sig","env_abs_matched_med")],
           setNames(HM[,c("cohort","sig","env_matched_med")], c("cohort","sig","env_abs_matched_med")))

lbl <- c(GGI="GGI\n(uniform)", MammaPrint="70-gene-like\n(near-uniform)",
         CorePAM="24-gene\n(concentrated)", OncotypeDX="21-gene-like\n(concentrated)")
pal <- c("#1b6ca8","#4ea1d3","#d95f02","#a33b1e"); names(pal) <- lbl

# --- Painel A: envelope com o painel PAREADO em 24 genes ---
md <- rbind(
  data.frame(sig=M$sig, env=M$env_abs_matched_med, cohort=M$cohort),
  data.frame(sig=A$sig[A$sig %in% c("CorePAM","OncotypeDX")],
             env=A$env_abs[A$sig %in% c("CorePAM","OncotypeDX")],
             cohort=A$cohort[A$sig %in% c("CorePAM","OncotypeDX")]))
md$lab <- factor(lbl[md$sig], levels=lbl)
pA <- ggplot(md, aes(lab, env, fill=lab)) +
  geom_boxplot(width=.55, alpha=.85, outlier.shape=NA) +
  geom_jitter(width=.10, size=1.3, alpha=.75) +
  scale_fill_manual(values=pal) +
  scale_y_continuous(breaks=seq(0,24,4)) +
  coord_cartesian(ylim=c(0,24)) +
  labs(title="A  Matched panel size (24 genes)",
       subtitle="Large panels randomly reduced to 24 genes;\nsmall scores at their native size",
       x=NULL, y="Genes removable before C-index drops 0.02  (of 24)") +
  theme_classic(base_size=10) +
  theme(legend.position="none", plot.subtitle=element_text(size=8, colour="grey30"),
        axis.text.x=element_text(size=7.2))

# --- Painel B: envelope absoluto x numero efetivo de genes (y LINEAR) ---
A$lab <- factor(lbl[A$sig], levels=lbl)
pB <- ggplot(A, aes(n_eff, env_abs, colour=lab, shape=lab)) +
  geom_point(size=2.4, alpha=.9) +
  scale_colour_manual(values=pal) +
  scale_x_continuous(breaks=c(10,20,40,60,80)) +
  scale_y_continuous(breaks=seq(0,90,15)) +
  labs(title="B  Envelope versus effective number of genes",
       subtitle="Each point is one signature in one cohort (six cohorts)",
       x=expression("Effective number of genes  "*n[eff]),
       y="Genes removable before C-index drops 0.02", colour=NULL, shape=NULL) +
  theme_classic(base_size=10) +
  theme(legend.position="right", legend.text=element_text(size=7.2),
        legend.key.height=unit(1.05,"lines"),
        plot.subtitle=element_text(size=8, colour="grey30"))

p <- pA + pB + plot_layout(widths=c(1,1.4))
ggsave(file.path(OUT,"Figure2_matched_neff.png"), p, width=10.2, height=4.2, dpi=300)
ggsave(file.path(OUT,"Figure2_matched_neff.pdf"), p, width=10.2, height=4.2)
ggsave(file.path(OUT,"Figure2_matched_neff.tiff"), p, width=10.2, height=4.2, dpi=300, compression="lzw")
cat("pontos no painel A:", nrow(md), " | painel B:", nrow(A), "\n")
cat("zeros no painel B:", sum(A$env_abs==0), "(todos devem aparecer)\n")
cat("DONE fig2\n")
