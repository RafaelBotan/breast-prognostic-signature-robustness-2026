# Part of the reproducibility repository for:
# "Robustness and transportability of breast-cancer prognostic gene signatures
#  under structured gene loss and platform change."
# Public cohorts only (GEO/TCGA/cBioPortal); no patient-level data are redistributed.
# Set COHORT_DIR to the folder with the processed cohort matrices (see README).

# LOCKED STRESS TEST v1 — degradacao sob perda estruturada de genes
# Travas: scoring rescalado PRE-ESPECIFICADO por assinatura (principio 7);
#         reportar distribuicao POR COORTE, sem numero unico (principio 6).
# Scoring uniforme = regra congelada do CorePAM: sum(z*w present)/sum(|w present|), high=pior.
suppressMessages({library(genefu); library(arrow); library(survival)})
set.seed(42)
base <- Sys.getenv("COHORT_DIR", "data/cohorts")
OUT  <- "results"
cohorts <- c("SCANB","TCGA_BRCA","METABRIC","GSE20685","GSE1456")
log <- function(...) cat(sprintf(...), "\n")

## ---- 1) Assinaturas: gene (HUGO) + peso, orientado high=pior ----
corepam <- c(EXO1=0.2174,NAT1=-0.2052,BLVRA=-0.1848,ACTR3B=-0.1405,MIA=-0.1193,
  MYBL2=-0.1183,PTTG1=0.1140,MDM2=-0.1074,SFRP1=-0.0907,GPR160=-0.0841,FOXC1=0.0717,
  PHGDH=0.0620,MYC=-0.0424,ESR1=0.0396,KRT5=-0.0277,CXXC5=0.0266,PGR=-0.0262,
  KRT17=-0.0163,CENPF=0.0150,FGFR4=0.0110,BCL2=-0.0105,MLPH=-0.0099,ERBB2=-0.0082,GRB7=-0.0053)

data(sig.oncotypedx); data(sig.ggi); data(sig.gene70)
# OncotypeDX: weight x sinal do grupo (estrogen = protetor; reference dropado)
od <- sig.oncotypedx
od <- od[od$group != "reference", ]
od_sign <- ifelse(od$group == "estrogen", -1, 1)
oncotype <- setNames(od$weight * od_sign, od$symbol)
# GGI: +1 se grade alto (3), -1 se baixo; usa HUGO, dropa NA
gg <- sig.ggi[!is.na(sig.ggi$HUGO.gene.symbol), ]
ggi <- setNames(ifelse(gg$grade == max(gg$grade,na.rm=TRUE), 1, -1), gg$HUGO.gene.symbol)
ggi <- ggi[!duplicated(names(ggi))]
# MammaPrint (gene70): peso = -correlacao (corr a bom prognostico -> high=pior)
g70 <- sig.gene70[!is.na(sig.gene70$HUGO.gene.symbol), ]
mammaprint <- setNames(-g70$correlation, g70$HUGO.gene.symbol)
mammaprint <- mammaprint[!duplicated(names(mammaprint))]

sigs <- list(CorePAM=corepam, OncotypeDX=oncotype, GGI=ggi, MammaPrint=mammaprint)
for(s in names(sigs)) log("  sig %s: %d genes", s, length(sigs[[s]]))

## ---- 2) Scoring congelado (principio 7) ----
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

## ---- 3) Loop coortes ----
rows <- list(); env <- list()
for(co in cohorts){
  ex <- as.data.frame(read_parquet(sprintf("%s/%s/expression_genelevel_preZ.parquet",base,co)))
  g <- ex$gene; ex$gene <- NULL; ex <- as.matrix(ex); rownames(ex) <- g
  ex <- ex[, !grepl("REPL$", colnames(ex)), drop=FALSE]
  zmat <- t(scale(t(ex)))                       # z-score por gene (linha)
  zmat <- zmat[is.finite(rowSums(zmat)), , drop=FALSE]
  cl <- as.data.frame(read_parquet(sprintf("%s/%s/clinical_FINAL.parquet",base,co)))
  key <- if(sum(colnames(zmat) %in% cl$patient_id) > sum(colnames(zmat) %in% cl$sample_id)) "patient_id" else "sample_id"
  m <- match(colnames(zmat), cl[[key]])
  tt <- cl$os_time_months[m]; ev <- cl$os_event[m]
  log("== %s == genes:%d amostras:%d (chave %s) eventos:%d", co, nrow(zmat), ncol(zmat), key, sum(ev,na.rm=TRUE))

  for(sg in names(sigs)){
    w0 <- sigs[[sg]]
    pres0 <- names(w0)[names(w0) %in% rownames(zmat)]
    if(length(pres0) < 5){ log("   %s: so %d/%d genes presentes -> pula", sg, length(pres0), length(w0)); next }
    w <- w0[pres0]
    # orientacao: full C-index; se <0.5 inverte sinal (mantem p/ todos os niveis)
    c_full <- cidx(tt,ev, score_sig(zmat,w))
    if(is.na(c_full)){ next }
    if(c_full < 0.5){ w <- -w; c_full <- cidx(tt,ev, score_sig(zmat,w)) }
    ng <- length(w)
    ord_w  <- names(sort(abs(w), decreasing=TRUE))            # remove maior |peso| 1o
    meanexpr <- rowMeans(ex[pres0,,drop=FALSE], na.rm=TRUE)
    ord_e  <- names(sort(meanexpr, decreasing=FALSE))          # remove menor expressao 1o
    # niveis de remocao
    levs <- 0:(ng-3)
    add <- function(lt, nrem, cval){
      rows[[length(rows)+1]] <<- data.frame(cohort=co, sig=sg, loss=lt,
        n_genes=ng, n_removed=nrem, frac_retained=round((ng-nrem)/ng,3),
        cindex=round(cval,4), drift=round(c_full-cval,4), stringsAsFactors=FALSE)
    }
    add("full",0,c_full)
    for(nrem in levs[-1]){
      # weight-directed
      keepw <- setdiff(names(w), ord_w[seq_len(nrem)])
      add("weight", nrem, cidx(tt,ev, score_sig(zmat, w[keepw])))
      # low-expression-directed
      keepe <- setdiff(names(w), ord_e[seq_len(nrem)])
      add("lowexpr", nrem, cidx(tt,ev, score_sig(zmat, w[keepe])))
      # random (B=20)
      cs <- replicate(20, { keepr <- sample(names(w), ng-nrem); cidx(tt,ev, score_sig(zmat, w[keepr])) })
      add("random", nrem, mean(cs, na.rm=TRUE))
    }
    # envelope: max fracao removida antes de drift >= 0.02 (por loss type)
    for(lt in c("weight","lowexpr","random")){
      sub <- do.call(rbind, rows); sub <- sub[sub$cohort==co & sub$sig==sg & sub$loss==lt,]
      bad <- sub[is.finite(sub$drift) & sub$drift>=0.02,]
      maxrem <- if(nrow(bad)==0) max(sub$n_removed) else min(bad$n_removed)-1
      env[[length(env)+1]] <- data.frame(cohort=co, sig=sg, loss=lt, n_genes=ng,
        c_full=round(c_full,4), max_removivel=maxrem, frac_tolerada=round(maxrem/ng,3),
        stringsAsFactors=FALSE)
    }
    log("   %s: full C=%.3f | env(weight)=%s genes", sg, c_full,
        env[[length(env)]]$max_removivel)
  }
}

curves <- do.call(rbind, rows); envelope <- do.call(rbind, env)
write.csv(curves,  file.path(OUT,"stress_curves.csv"),  row.names=FALSE)
write.csv(envelope,file.path(OUT,"stress_envelope.csv"),row.names=FALSE)
log("\nsalvos stress_curves.csv (%d linhas) e stress_envelope.csv (%d linhas)", nrow(curves), nrow(envelope))

## ---- 4) Figura: degradacao weight-directed por coorte ----
ok <- requireNamespace("ggplot2", quietly=TRUE)
if(ok){
  library(ggplot2)
  d <- curves[curves$loss %in% c("weight","full"),]
  p <- ggplot(d, aes(frac_retained, cindex, color=sig)) +
    geom_line(linewidth=0.8) + geom_point(size=0.7) +
    facet_wrap(~cohort, scales="free_y") +
    scale_x_reverse() +
    labs(title="Locked stress test: degradacao por perda dirigida por peso",
         subtitle="C-index (Harrell) vs fracao de genes retidos | por coorte (principio 6: distribuicao, nao numero unico)",
         x="fracao de genes retidos (->menos genes)", y="C-index", color="assinatura") +
    theme_classic(base_size=11)
  ggsave(file.path(OUT,"stress_degradation.png"), p, width=11, height=7, dpi=130)
  log("salvo stress_degradation.png")
}
cat("DONE\n")
