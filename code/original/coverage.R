# Part of the reproducibility repository for:
# "Robustness and transportability of breast-cancer prognostic gene signatures
#  under structured gene loss and platform change."
# Public cohorts only (GEO/TCGA/cBioPortal); no patient-level data are redistributed.
# Set COHORT_DIR to the folder with the processed cohort matrices (see README).

suppressMessages({library(genefu); library(arrow)})
data(sig.ggi); data(sig.gene70); data(sig.oncotypedx); data(sig.gene76); data(pam50)
base <- Sys.getenv("COHORT_DIR", "data/cohorts")
cohorts <- c("SCANB","TCGA_BRCA","METABRIC","GSE20685","GSE1456")
corepam24 <- c("EXO1","PTTG1","CENPF","MYBL2","MYC","FOXC1","KRT5","KRT17","ERBB2","GRB7","MIA","ESR1","PGR","BCL2","NAT1","MLPH","BLVRA","ACTR3B","MDM2","SFRP1","GPR160","PHGDH","CXXC5","FGFR4")
siggenes <- list(
  OncotypeDX = unique(sig.oncotypedx$symbol),
  CorePAM    = corepam24,
  ROR_PAM50  = rownames(pam50$centroids),
  MammaPrint = unique(sig.gene70$NCBI.gene.symbol),
  gene76     = unique(sig.gene76$NCBI.gene.symbol),
  GGI        = unique(sig.ggi$NCBI.gene.symbol)
)
# harmonizar 1 sinonimo conhecido p/ ser justo
harmonize_symbol <- function(v){ v[v=="CTSL2"]<-"CTSV"; v }
genes_co <- lapply(cohorts, function(co){
  g <- as.data.frame(read_parquet(sprintf("%s/%s/expression_genelevel_preZ.parquet",base,co)))$gene
  g })
names(genes_co) <- cohorts
cov <- sapply(cohorts, function(co){
  G <- genes_co[[co]]
  sapply(names(siggenes), function(s){
    sg <- harmonize_symbol(siggenes[[s]]); sg<-sg[!is.na(sg)&sg!=""]
    round(100*mean(sg %in% G),0)
  })
})
cat("=== COBERTURA (% genes da assinatura presentes na coorte) ===\n")
cat(sprintf("%-11s n_genes ", "")); cat(paste(sprintf("%9s",cohorts),collapse=""),"\n")
for(s in names(siggenes)){
  cat(sprintf("%-11s %5d  ", s, length(unique(harmonize_symbol(siggenes[[s]])))))
  cat(paste(sprintf("%8d%%", cov[s,]),collapse=""),"\n")
}
write.csv(cov, "results/coverage.csv")
