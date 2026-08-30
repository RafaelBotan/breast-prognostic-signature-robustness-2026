# STRESS TEST v2 — responde a objecao da conselheira (2026-07-29):
#   "a relacao redundancia->robustez pode ser MECANICA porque o envelope
#    e expresso como FRACAO de genes removidos; painel maior tolera fracao
#    maior por construcao. Testar tamanho pareado / numero EFETIVO de genes
#    e a incerteza das curvas."
#
# Tres experimentos:
#   A) ENVELOPE EM GENES ABSOLUTOS (ja calculavel do v1, refeito aqui p/ completude)
#   B) TAMANHO PAREADO: subamostrar GGI e MammaPrint para 24 genes (= tamanho do
#      CorePAM) e remedir o envelope. Se o painel grande subamostrado CONTINUA
#      robusto -> o motor e a DISTRIBUICAO DE PESO, nao o tamanho.
#   C) NUMERO EFETIVO DE GENES (participation ratio) vs envelope, + IC bootstrap
#      do envelope (reamostragem de PACIENTES).
suppressMessages({library(genefu); library(arrow); library(survival)})
set.seed(20260729)
# Portable path adaptation for the public repository; analytic logic is unchanged.
base <- Sys.getenv("COHORT_DIR", "data/cohorts")
OUT  <- Sys.getenv("REVISION_OUTPUT_DIR", "results")
cohorts <- c("SCANB","TCGA_BRCA","METABRIC","GSE20685","GSE1456")
log <- function(...) { cat(sprintf(...), "\n"); flush.console() }

## ---- assinaturas (identicas ao v1) ----
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

## ---- metricas de concentracao de peso ----
gini <- function(w){ a <- sort(abs(w)); n <- length(a)
  if(sum(a)==0) return(NA_real_); sum((2*seq_len(n)-n-1)*a)/(n*sum(a)) }
# numero EFETIVO de genes (participation ratio / inverso de Simpson sobre peso^2)
neff <- function(w){ a <- abs(w); if(sum(a)==0) return(NA_real_); (sum(a)^2)/sum(a^2) }

score_sig <- function(zmat, w){
  pres <- names(w)[names(w) %in% rownames(zmat)]
  if(length(pres) < 3) return(rep(NA_real_, ncol(zmat)))
  wp <- w[pres]
  as.numeric(colSums(zmat[pres,,drop=FALSE] * wp) / sum(abs(wp)))
}
cidx <- function(time,event,sc){
  ok <- is.finite(time)&is.finite(event)&is.finite(sc)
  if(sum(ok)<20 || length(unique(sc[ok]))<3) return(NA_real_)
  tryCatch(concordance(Surv(time[ok],event[ok]) ~ sc[ok])$concordance, error=function(e) NA_real_)
}
# envelope com SAIDA ANTECIPADA: quantos genes da p/ remover (weight-directed)
# antes do C-index cair >= thr. Retorna n absoluto de genes removiveis.
env_weight <- function(zmat,tt,ev,w,thr=0.02){
  c_full <- cidx(tt,ev,score_sig(zmat,w))
  if(!is.finite(c_full)) return(c(NA,NA))
  if(c_full < 0.5){ w <- -w; c_full <- cidx(tt,ev,score_sig(zmat,w)) }
  ng <- length(w); ord <- names(sort(abs(w), decreasing=TRUE))
  maxrem <- ng-3
  for(nrem in seq_len(ng-3)){
    keep <- setdiff(names(w), ord[seq_len(nrem)])
    cv <- cidx(tt,ev,score_sig(zmat,w[keep]))
    if(is.finite(cv) && (c_full-cv) >= thr){ maxrem <- nrem-1; break }
  }
  c(c_full, maxrem)
}

B_MATCH <- 30    # subconjuntos aleatorios no experimento de tamanho pareado
B_BOOT  <- 100   # replicas bootstrap de pacientes p/ IC do envelope
TARGET  <- 24    # tamanho pareado (= CorePAM)

resA <- list(); resB <- list(); resC <- list()

for(co in cohorts){
  ex <- as.data.frame(read_parquet(sprintf("%s/%s/expression_genelevel_preZ.parquet",base,co)))
  g <- ex$gene; ex$gene <- NULL; ex <- as.matrix(ex); rownames(ex) <- g
  ex <- ex[, !grepl("REPL$", colnames(ex)), drop=FALSE]
  zmat <- t(scale(t(ex))); zmat <- zmat[is.finite(rowSums(zmat)), , drop=FALSE]
  cl <- as.data.frame(read_parquet(sprintf("%s/%s/clinical_FINAL.parquet",base,co)))
  key <- if(sum(colnames(zmat) %in% cl$patient_id) > sum(colnames(zmat) %in% cl$sample_id)) "patient_id" else "sample_id"
  m <- match(colnames(zmat), cl[[key]])
  tt <- cl$os_time_months[m]; ev <- cl$os_event[m]
  log("== %s == genes:%d amostras:%d eventos:%d", co, nrow(zmat), ncol(zmat), sum(ev,na.rm=TRUE))

  for(sg in names(sigs)){
    w0 <- sigs[[sg]]; pres0 <- names(w0)[names(w0) %in% rownames(zmat)]
    if(length(pres0) < 5) next
    w <- w0[pres0]; ng <- length(w)
    e <- env_weight(zmat,tt,ev,w); c_full <- e[1]; maxrem <- e[2]
    if(!is.finite(c_full)) next
    # ---- A + C: envelope absoluto, fracao, concentracao de peso ----
    resA[[length(resA)+1]] <- data.frame(cohort=co, sig=sg, n_genes=ng,
      n_genes_nominal=length(w0), c_full=round(c_full,4),
      env_abs=maxrem, env_frac=round(maxrem/ng,4),
      gini=round(gini(w),4), n_eff=round(neff(w),2),
      stringsAsFactors=FALSE)
    log("   %s: n=%d n_eff=%.1f gini=%.2f | C=%.3f env_abs=%d (%.1f%%)",
        sg, ng, neff(w), gini(w), c_full, maxrem, 100*maxrem/ng)

    # ---- C: IC bootstrap do envelope (reamostra PACIENTES) ----
    nsamp <- ncol(zmat); bs <- integer(0)
    for(b in seq_len(B_BOOT)){
      idx <- sample.int(nsamp, nsamp, replace=TRUE)
      eb <- env_weight(zmat[,idx,drop=FALSE], tt[idx], ev[idx], w)
      if(is.finite(eb[2])) bs <- c(bs, eb[2])
    }
    if(length(bs) >= 20){
      q <- quantile(bs, c(0.025,0.975), na.rm=TRUE)
      resC[[length(resC)+1]] <- data.frame(cohort=co, sig=sg, n_genes=ng,
        env_abs=maxrem, boot_med=median(bs), boot_lo=as.numeric(q[1]), boot_hi=as.numeric(q[2]),
        boot_lo_frac=round(as.numeric(q[1])/ng,4), boot_hi_frac=round(as.numeric(q[2])/ng,4),
        n_boot=length(bs), stringsAsFactors=FALSE)
      log("      boot env_abs mediana=%.0f IC95%% [%.0f, %.0f]", median(bs), q[1], q[2])
    }

    # ---- B: TAMANHO PAREADO (so p/ paineis maiores que o alvo) ----
    if(ng > TARGET + 4){
      envs <- numeric(0); cfs <- numeric(0)
      for(b in seq_len(B_MATCH)){
        sub <- sample(names(w), TARGET)
        eb <- env_weight(zmat,tt,ev,w[sub])
        if(is.finite(eb[2])){ envs <- c(envs, eb[2]); cfs <- c(cfs, eb[1]) }
      }
      if(length(envs) >= 10){
        resB[[length(resB)+1]] <- data.frame(cohort=co, sig=sg, n_full=ng, n_matched=TARGET,
          env_abs_full=maxrem, env_abs_matched_med=median(envs),
          env_abs_matched_lo=as.numeric(quantile(envs,0.1)),
          env_abs_matched_hi=as.numeric(quantile(envs,0.9)),
          env_frac_matched_med=round(median(envs)/TARGET,4),
          c_full_matched_med=round(median(cfs),4), n_draws=length(envs),
          stringsAsFactors=FALSE)
        log("      PAREADO(24): env_abs mediana=%.1f (P10-P90 %.0f-%.0f) | frac=%.1f%% | C=%.3f",
            median(envs), quantile(envs,0.1), quantile(envs,0.9), 100*median(envs)/TARGET, median(cfs))
      }
    }
  }
}

A <- do.call(rbind,resA); B <- do.call(rbind,resB); C <- do.call(rbind,resC)
write.csv(A, file.path(OUT,"v2_envelope_abs_neff.csv"), row.names=FALSE)
write.csv(B, file.path(OUT,"v2_matched_size.csv"),     row.names=FALSE)
write.csv(C, file.path(OUT,"v2_bootstrap_env.csv"),    row.names=FALSE)
log("\n== salvos v2_envelope_abs_neff.csv (%d) / v2_matched_size.csv (%d) / v2_bootstrap_env.csv (%d)",
    nrow(A), nrow(B), nrow(C))

## ---- relacao envelope ~ numero efetivo de genes ----
if(!is.null(A) && nrow(A) > 5){
  ct <- suppressWarnings(cor.test(A$n_eff, A$env_abs, method="spearman"))
  ct2<- suppressWarnings(cor.test(A$n_genes, A$env_abs, method="spearman"))
  log("Spearman env_abs ~ n_eff   : rho=%.3f p=%.2g", ct$estimate, ct$p.value)
  log("Spearman env_abs ~ n_genes : rho=%.3f p=%.2g", ct2$estimate, ct2$p.value)
}
cat("DONE v2\n")
