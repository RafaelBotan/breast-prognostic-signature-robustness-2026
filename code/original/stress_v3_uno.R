# STRESS TEST v3 — sensibilidade do envelope ao estimador de concordancia.
# O v1/v2 usaram o C de Harrell, que e sensivel ao padrao de censura.
# Aqui repetimos o envelope com o C de UNO (ponderado por 1/G(t)^2, robusto a
# censura dependente do tempo), disponivel de graca em survival::concordance
# via timewt="n/G2". Se a ordenacao de robustez se mantiver, cai a limitacao
# "Uno nao foi computado" que constava do manuscrito v1.
suppressMessages({library(genefu); library(arrow); library(survival)})
set.seed(20260729)
# Portable path adaptation for the public repository; analytic logic is unchanged.
base <- Sys.getenv("COHORT_DIR", "data/cohorts")
OUT  <- Sys.getenv("REVISION_OUTPUT_DIR", "results")
cohorts <- c("SCANB","TCGA_BRCA","METABRIC","GSE20685","GSE1456")
log <- function(...) { cat(sprintf(...), "\n"); flush.console() }

corepam <- c(EXO1=0.2174,NAT1=-0.2052,BLVRA=-0.1848,ACTR3B=-0.1405,MIA=-0.1193,
  MYBL2=-0.1183,PTTG1=0.1140,MDM2=-0.1074,SFRP1=-0.0907,GPR160=-0.0841,FOXC1=0.0717,
  PHGDH=0.0620,MYC=-0.0424,ESR1=0.0396,KRT5=-0.0277,CXXC5=0.0266,PGR=-0.0262,
  KRT17=-0.0163,CENPF=0.0150,FGFR4=0.0110,BCL2=-0.0105,MLPH=-0.0099,ERBB2=-0.0082,GRB7=-0.0053)
data(sig.oncotypedx); data(sig.ggi); data(sig.gene70)
od <- sig.oncotypedx; od <- od[od$group != "reference", ]
oncotype <- setNames(od$weight * ifelse(od$group=="estrogen",-1,1), od$symbol)
gg  <- sig.ggi[!is.na(sig.ggi$HUGO.gene.symbol), ]
ggi <- setNames(ifelse(gg$grade==max(gg$grade,na.rm=TRUE),1,-1), gg$HUGO.gene.symbol)
ggi <- ggi[!duplicated(names(ggi))]
g70 <- sig.gene70[!is.na(sig.gene70$HUGO.gene.symbol), ]
mammaprint <- setNames(-g70$correlation, g70$HUGO.gene.symbol)
mammaprint <- mammaprint[!duplicated(names(mammaprint))]
sigs <- list(CorePAM=corepam, OncotypeDX=oncotype, GGI=ggi, MammaPrint=mammaprint)

score_sig <- function(zmat, w){
  pres <- names(w)[names(w) %in% rownames(zmat)]
  if(length(pres) < 3) return(rep(NA_real_, ncol(zmat)))
  wp <- w[pres]; as.numeric(colSums(zmat[pres,,drop=FALSE] * wp) / sum(abs(wp)))
}
# estimador parametrizavel: Harrell (default) ou Uno (timewt="n/G2")
cidx <- function(time,event,sc,wt=c("n","n/G2")){
  wt <- match.arg(wt)
  ok <- is.finite(time)&is.finite(event)&is.finite(sc)
  if(sum(ok)<20 || length(unique(sc[ok]))<3) return(NA_real_)
  tryCatch(concordance(Surv(time[ok],event[ok]) ~ sc[ok], timewt=wt)$concordance,
           error=function(e) NA_real_)
}
env_weight <- function(zmat,tt,ev,w,thr=0.02,wt="n"){
  c_full <- cidx(tt,ev,score_sig(zmat,w),wt)
  if(!is.finite(c_full)) return(c(NA,NA))
  if(c_full < 0.5){ w <- -w; c_full <- cidx(tt,ev,score_sig(zmat,w),wt) }
  ng <- length(w); ord <- names(sort(abs(w), decreasing=TRUE)); maxrem <- ng-3
  for(nrem in seq_len(ng-3)){
    keep <- setdiff(names(w), ord[seq_len(nrem)])
    cv <- cidx(tt,ev,score_sig(zmat,w[keep]),wt)
    if(is.finite(cv) && (c_full-cv) >= thr){ maxrem <- nrem-1; break }
  }
  c(c_full, maxrem)
}

TARGET <- 24; B_MATCH <- 30
res <- list(); resM <- list()
for(co in cohorts){
  ex <- as.data.frame(read_parquet(sprintf("%s/%s/expression_genelevel_preZ.parquet",base,co)))
  g <- ex$gene; ex$gene <- NULL; ex <- as.matrix(ex); rownames(ex) <- g
  ex <- ex[, !grepl("REPL$", colnames(ex)), drop=FALSE]
  zmat <- t(scale(t(ex))); zmat <- zmat[is.finite(rowSums(zmat)), , drop=FALSE]
  cl <- as.data.frame(read_parquet(sprintf("%s/%s/clinical_FINAL.parquet",base,co)))
  key <- if(sum(colnames(zmat) %in% cl$patient_id) > sum(colnames(zmat) %in% cl$sample_id)) "patient_id" else "sample_id"
  m <- match(colnames(zmat), cl[[key]]); tt <- cl$os_time_months[m]; ev <- cl$os_event[m]
  log("== %s ==", co)
  for(sg in names(sigs)){
    w0 <- sigs[[sg]]; pres0 <- names(w0)[names(w0) %in% rownames(zmat)]
    if(length(pres0) < 5) next
    w <- w0[pres0]; ng <- length(w)
    eH <- env_weight(zmat,tt,ev,w,wt="n")
    eU <- env_weight(zmat,tt,ev,w,wt="n/G2")
    if(!is.finite(eH[1])) next
    res[[length(res)+1]] <- data.frame(cohort=co, sig=sg, n_genes=ng,
      c_harrell=round(eH[1],4), env_harrell=eH[2], env_harrell_frac=round(eH[2]/ng,4),
      c_uno=round(eU[1],4), env_uno=eU[2], env_uno_frac=round(eU[2]/ng,4),
      stringsAsFactors=FALSE)
    log("   %-11s n=%2d | Harrell C=%.3f env=%2d (%.0f%%) | Uno C=%.3f env=%2d (%.0f%%)",
        sg, ng, eH[1], eH[2], 100*eH[2]/ng, eU[1], eU[2], 100*eU[2]/ng)
    if(ng > TARGET + 4){
      envs <- numeric(0)
      for(b in seq_len(B_MATCH)){
        sub <- sample(names(w), TARGET)
        eb <- env_weight(zmat,tt,ev,w[sub],wt="n/G2")
        if(is.finite(eb[2])) envs <- c(envs, eb[2])
      }
      if(length(envs) >= 10){
        resM[[length(resM)+1]] <- data.frame(cohort=co, sig=sg, n_full=ng,
          env_uno_matched_med=median(envs), env_uno_matched_lo=as.numeric(quantile(envs,0.1)),
          env_uno_matched_hi=as.numeric(quantile(envs,0.9)), n_draws=length(envs),
          stringsAsFactors=FALSE)
        log("      PAREADO(24) sob Uno: env mediana=%.0f (P10-P90 %.0f-%.0f)",
            median(envs), quantile(envs,0.1), quantile(envs,0.9))
      }
    }
  }
}
R <- do.call(rbind,res); M <- do.call(rbind,resM)
write.csv(R, file.path(OUT,"v3_uno_vs_harrell.csv"), row.names=FALSE)
write.csv(M, file.path(OUT,"v3_uno_matched.csv"),    row.names=FALSE)
log("\nsalvos v3_uno_vs_harrell.csv (%d) e v3_uno_matched.csv (%d)", nrow(R), nrow(M))
if(nrow(R)>5){
  ct <- suppressWarnings(cor.test(R$env_harrell, R$env_uno, method="spearman"))
  log("Concordancia da ORDENACAO Harrell vs Uno: Spearman rho=%.3f p=%.2g", ct$estimate, ct$p.value)
}
cat("DONE v3\n")
