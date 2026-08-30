# STRESS TEST v4 — (a) repete na coorte HELD-OUT GSE7390 as analises novas do v2/v3
# (tamanho pareado, numero efetivo de genes, IC bootstrap, C de Uno) e
# (b) gera a Figura 2 do manuscrito reenquadrado.
suppressMessages({library(GEOquery); library(genefu); library(survival)})
set.seed(20260729)
# Portable path adaptation for the public repository; analytic logic is unchanged.
OUT <- Sys.getenv("REVISION_OUTPUT_DIR", "results")
log <- function(...) { cat(sprintf(...), "\n"); flush.console() }
options(timeout=900)

## ---------- assinaturas (identicas ao pipeline travado) ----------
corepam <- c(EXO1=0.2174,NAT1=-0.2052,BLVRA=-0.1848,ACTR3B=-0.1405,MIA=-0.1193,MYBL2=-0.1183,PTTG1=0.1140,
  MDM2=-0.1074,SFRP1=-0.0907,GPR160=-0.0841,FOXC1=0.0717,PHGDH=0.0620,MYC=-0.0424,ESR1=0.0396,KRT5=-0.0277,
  CXXC5=0.0266,PGR=-0.0262,KRT17=-0.0163,CENPF=0.0150,FGFR4=0.0110,BCL2=-0.0105,MLPH=-0.0099,ERBB2=-0.0082,GRB7=-0.0053)
data(sig.oncotypedx); data(sig.ggi); data(sig.gene70)
od <- sig.oncotypedx; od <- od[od$group!="reference",]
oncotype <- setNames(od$weight*ifelse(od$group=="estrogen",-1,1), od$symbol)
gg <- sig.ggi[!is.na(sig.ggi$HUGO.gene.symbol),]
ggi <- setNames(ifelse(gg$grade==max(gg$grade,na.rm=TRUE),1,-1), gg$HUGO.gene.symbol); ggi <- ggi[!duplicated(names(ggi))]
g70 <- sig.gene70[!is.na(sig.gene70$HUGO.gene.symbol),]
mammaprint <- setNames(-g70$correlation, g70$HUGO.gene.symbol); mammaprint <- mammaprint[!duplicated(names(mammaprint))]
sigs <- list(CorePAM=corepam, OncotypeDX=oncotype, GGI=ggi, MammaPrint=mammaprint)

gini <- function(w){ a<-sort(abs(w)); n<-length(a); if(sum(a)==0) return(NA_real_); sum((2*seq_len(n)-n-1)*a)/(n*sum(a)) }
neff <- function(w){ a<-abs(w); if(sum(a)==0) return(NA_real_); (sum(a)^2)/sum(a^2) }
score_sig <- function(z,w){ p<-names(w)[names(w)%in%rownames(z)]; if(length(p)<3) return(rep(NA_real_,ncol(z)))
  as.numeric(colSums(z[p,,drop=FALSE]*w[p])/sum(abs(w[p]))) }
cidx <- function(t,e,s,wt="n"){ ok<-is.finite(t)&is.finite(e)&is.finite(s)
  if(sum(ok)<20||length(unique(s[ok]))<3) return(NA_real_)
  tryCatch(concordance(Surv(t[ok],e[ok])~s[ok], timewt=wt)$concordance, error=function(x) NA_real_) }
env_weight <- function(z,tt,ev,w,thr=0.02,wt="n"){
  cf <- cidx(tt,ev,score_sig(z,w),wt); if(!is.finite(cf)) return(c(NA,NA))
  if(cf<0.5){ w <- -w; cf <- cidx(tt,ev,score_sig(z,w),wt) }
  ng<-length(w); ord<-names(sort(abs(w),decreasing=TRUE)); mx<-ng-3
  for(nr in seq_len(ng-3)){
    cv <- cidx(tt,ev,score_sig(z,w[setdiff(names(w),ord[seq_len(nr)])]),wt)
    if(is.finite(cv) && (cf-cv)>=thr){ mx<-nr-1; break }
  }
  c(cf,mx)
}

## ---------- (a) HELD-OUT GSE7390 ----------
log("baixando GSE7390...")
gse <- tryCatch(getGEO("GSE7390", GSEMatrix=TRUE, AnnotGPL=TRUE, destdir=tempdir())[[1]],
                error=function(e){ log("ERRO: %s", conditionMessage(e)); NULL })
heldout_ok <- FALSE
if(!is.null(gse)){
  ex<-exprs(gse); fd<-fData(gse); pd<-pData(gse)
  sym_col <- colnames(fd)[grep("symbol",colnames(fd),ignore.case=TRUE)][1]
  findcol <- function(pat){ c<-colnames(pd)[grep(pat,colnames(pd),ignore.case=TRUE)]; if(length(c)) c[1] else NA }
  tcol<-findcol("t.dmfs|t_dmfs|dmfs.time|time.dmfs"); ecol<-findcol("e.dmfs|e_dmfs|dmfs.event|event.dmfs")
  pnum <- function(x) as.numeric(gsub(".*: *","",as.character(x)))
  if(is.na(tcol)||is.na(ecol)){
    ch <- pd[,grep("characteristics",colnames(pd)),drop=FALSE]; flat<-apply(ch,2,as.character)
    t2<-which(apply(flat,2,function(c) any(grepl("t.dmfs|tdm",c,ignore.case=TRUE))))
    e2<-which(apply(flat,2,function(c) any(grepl("e.dmfs|edm",c,ignore.case=TRUE))))
    time<-pnum(ch[,t2[1]]); event<-pnum(ch[,e2[1]])
  } else { time<-pnum(pd[[tcol]]); event<-pnum(pd[[ecol]]) }
  sym<-sub(" ///.*","",fd[[sym_col]]); keep<-sym!=""&!is.na(sym); ex<-ex[keep,]; sym<-sym[keep]
  o<-order(apply(ex,1,var,na.rm=TRUE),decreasing=TRUE); ex<-ex[o,]; sym<-sym[o]
  dup<-duplicated(sym); exg<-ex[!dup,]; rownames(exg)<-sym[!dup]
  z<-t(scale(t(exg))); z<-z[is.finite(rowSums(z)),,drop=FALSE]
  log("GSE7390: %d genes, %d amostras, %d eventos", nrow(z), ncol(z), sum(event==1,na.rm=TRUE))
  rows<-list(); rowsM<-list()
  for(sg in names(sigs)){
    w0<-sigs[[sg]]; pres<-names(w0)[names(w0)%in%rownames(z)]
    if(length(pres)<5) next
    w<-w0[pres]; ng<-length(w)
    eH<-env_weight(z,time,event,w,wt="n"); eU<-env_weight(z,time,event,w,wt="n/G2")
    if(!is.finite(eH[1])) next
    bs<-integer(0); ns<-ncol(z)
    for(b in 1:100){ i<-sample.int(ns,ns,replace=TRUE)
      eb<-env_weight(z[,i,drop=FALSE],time[i],event[i],w); if(is.finite(eb[2])) bs<-c(bs,eb[2]) }
    q<-if(length(bs)>=20) quantile(bs,c(.025,.975),na.rm=TRUE) else c(NA,NA)
    rows[[length(rows)+1]]<-data.frame(cohort="GSE7390",sig=sg,n_genes=ng,n_eff=round(neff(w),2),
      gini=round(gini(w),4),c_harrell=round(eH[1],4),env_abs=eH[2],env_frac=round(eH[2]/ng,4),
      boot_lo=as.numeric(q[1]),boot_hi=as.numeric(q[2]),c_uno=round(eU[1],4),env_uno=eU[2],
      stringsAsFactors=FALSE)
    log("  %-11s n=%2d n_eff=%.1f gini=%.2f | C=%.3f env=%2d (%.0f%%) IC[%.0f,%.0f] | Uno C=%.3f env=%d",
        sg,ng,neff(w),gini(w),eH[1],eH[2],100*eH[2]/ng,q[1],q[2],eU[1],eU[2])
    if(ng > 28){
      envs<-numeric(0)
      for(b in 1:30){ sub<-sample(names(w),24); eb<-env_weight(z,time,event,w[sub]); if(is.finite(eb[2])) envs<-c(envs,eb[2]) }
      if(length(envs)>=10){
        rowsM[[length(rowsM)+1]]<-data.frame(cohort="GSE7390",sig=sg,n_full=ng,
          env_matched_med=median(envs),env_matched_lo=as.numeric(quantile(envs,.1)),
          env_matched_hi=as.numeric(quantile(envs,.9)),n_draws=length(envs),stringsAsFactors=FALSE)
        log("     PAREADO(24): mediana=%.0f (P10-P90 %.0f-%.0f)",median(envs),quantile(envs,.1),quantile(envs,.9))
      }
    }
  }
  if(length(rows)){ write.csv(do.call(rbind,rows), file.path(OUT,"v4_heldout_GSE7390.csv"), row.names=FALSE); heldout_ok<-TRUE }
  if(length(rowsM)) write.csv(do.call(rbind,rowsM), file.path(OUT,"v4_heldout_matched.csv"), row.names=FALSE)
}

## ---------- (b) FIGURA 2 ----------
ok <- requireNamespace("ggplot2",quietly=TRUE) && requireNamespace("patchwork",quietly=TRUE)
A <- read.csv(file.path(OUT,"v2_envelope_abs_neff.csv"), stringsAsFactors=FALSE)
M <- read.csv(file.path(OUT,"v2_matched_size.csv"),      stringsAsFactors=FALSE)
if(heldout_ok){
  H <- read.csv(file.path(OUT,"v4_heldout_GSE7390.csv"), stringsAsFactors=FALSE)
  A <- rbind(A[,c("cohort","sig","n_genes","env_abs","env_frac","gini","n_eff")],
             data.frame(cohort=H$cohort,sig=H$sig,n_genes=H$n_genes,env_abs=H$env_abs,
                        env_frac=H$env_frac,gini=H$gini,n_eff=H$n_eff))
}
lbl <- c(GGI="GGI (uniform weights)", MammaPrint="70-gene-like (near-uniform)",
         CorePAM="24-gene (concentrated)", OncotypeDX="21-gene-like (concentrated)")
A$lab <- lbl[A$sig]; A$lab <- factor(A$lab, levels=lbl)

if(ok){
  library(ggplot2); library(patchwork)
  # Painel A: envelope pareado a 24 genes
  md <- data.frame(
    cohort = c(M$cohort, A$cohort[A$sig %in% c("CorePAM")]),
    sig    = c(M$sig,    A$sig[A$sig %in% c("CorePAM")]),
    env    = c(M$env_abs_matched_med, A$env_abs[A$sig %in% c("CorePAM")]),
    stringsAsFactors = FALSE)
  md$lab <- lbl[md$sig]; md$lab <- factor(md$lab, levels=lbl)
  pA <- ggplot(md, aes(lab, env, fill=lab)) +
    geom_boxplot(outlier.size=0.6, width=0.55, alpha=0.85) +
    geom_jitter(width=0.09, size=1.1, alpha=0.7) +
    scale_y_continuous(limits=c(0,24), breaks=seq(0,24,4)) +
    labs(title="A. Matched panel size (24 genes)",
         subtitle="Genes removable before the C-index falls by 0.02",
         x=NULL, y="Genes removable (of 24)") +
    theme_classic(base_size=10) + theme(legend.position="none",
      axis.text.x=element_text(angle=18, hjust=1))
  # Painel B: envelope absoluto vs numero efetivo de genes
  pB <- ggplot(A, aes(n_eff, env_abs, color=lab, shape=lab)) +
    geom_point(size=2.2, alpha=0.85) +
    scale_x_log10() + scale_y_log10() +
    labs(title="B. Envelope versus effective number of genes",
         subtitle="Each point is one signature in one cohort",
         x=expression("Effective number of genes  "*n[eff]*"  (log scale)"),
         y="Genes removable (log scale)", color=NULL, shape=NULL) +
    theme_classic(base_size=10) + theme(legend.position="right",
      legend.text=element_text(size=7.5))
  p <- pA + pB + plot_layout(widths=c(1,1.35))
  ggsave(file.path(OUT,"Figure2_matched_neff.png"), p, width=10.5, height=4.3, dpi=300)
  ggsave(file.path(OUT,"Figure2_matched_neff.pdf"), p, width=10.5, height=4.3)
  log("salvo Figure2_matched_neff.png/.pdf")
} else log("ggplot2/patchwork ausente — figura 2 NAO gerada")
cat("DONE v4\n")
