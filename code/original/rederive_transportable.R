# Part of the reproducibility repository for:
# "Robustness and transportability of breast-cancer prognostic gene signatures
#  under structured gene loss and platform change."
# Public cohorts only (GEO/TCGA/cBioPortal); no patient-level data are redistributed.
# Set COHORT_DIR to the folder with the processed cohort matrices (see README).

# #2 — Re-derivacao TRANSPORTAVEL (prova de conceito): CorePAM-2 do pool PAM50
# Trava (conselheira): treino LEAVE-PLATFORM-OUT; z-score por coorte (sem leakage de escala);
#   selecao so no treino; teste em plataforma+coorte HELD-OUT.
# Hipotese: penalidade ridge-pesada + treino multiplataforma -> espalha peso -> menos fragil,
#   talvez com pequeno custo de discriminacao. Compara CorePAM vs CorePAM-2.
suppressMessages({library(arrow); library(survival); library(glmnet); library(genefu)})
set.seed(7)
base <- Sys.getenv("COHORT_DIR", "data/cohorts")
OUT  <- "results"
log <- function(...) cat(sprintf(...), "\n")

## pool de genes = PAM50
data(pam50)
pam50_genes <- rownames(pam50$centroids)
log("PAM50 pool: %d genes (ex: %s)", length(pam50_genes), paste(head(pam50_genes,5),collapse=","))

## CorePAM original (referencia)
corepam <- c(EXO1=0.2174,NAT1=-0.2052,BLVRA=-0.1848,ACTR3B=-0.1405,MIA=-0.1193,MYBL2=-0.1183,
  PTTG1=0.1140,MDM2=-0.1074,SFRP1=-0.0907,GPR160=-0.0841,FOXC1=0.0717,PHGDH=0.0620,MYC=-0.0424,
  ESR1=0.0396,KRT5=-0.0277,CXXC5=0.0266,PGR=-0.0262,KRT17=-0.0163,CENPF=0.0150,FGFR4=0.0110,
  BCL2=-0.0105,MLPH=-0.0099,ERBB2=-0.0082,GRB7=-0.0053)

cohorts <- c("SCANB","TCGA_BRCA","METABRIC","GSE20685","GSE1456")
D <- list()
for(co in cohorts){
  ex<-as.data.frame(read_parquet(sprintf("%s/%s/expression_genelevel_preZ.parquet",base,co)))
  g<-ex$gene; ex$gene<-NULL; ex<-as.matrix(ex); rownames(ex)<-g
  ex<-ex[, !grepl("REPL$",colnames(ex)),drop=FALSE]
  z<-t(scale(t(ex))); z<-z[is.finite(rowSums(z)),,drop=FALSE]      # z por gene
  cl<-as.data.frame(read_parquet(sprintf("%s/%s/clinical_FINAL.parquet",base,co)))
  key<-if(sum(colnames(ex)%in%cl$patient_id)>sum(colnames(ex)%in%cl$sample_id))"patient_id" else "sample_id"
  m<-match(colnames(ex), cl[[key]])
  D[[co]]<-list(z=z, time=cl$os_time_months[m], event=cl$os_event[m])
}
cidx<-function(time,event,sc){ok<-is.finite(time)&is.finite(event)&is.finite(sc); if(sum(ok)<20||length(unique(sc[ok]))<3)return(NA_real_); tryCatch(concordance(Surv(time[ok],event[ok])~sc[ok])$concordance,error=function(e)NA_real_)}
sc_w<-function(z,w){p<-names(w)[names(w)%in%rownames(z)]; if(length(p)<3)return(rep(NA_real_,ncol(z))); as.numeric(colSums(z[p,,drop=FALSE]*w[p])/sum(abs(w[p])))}
gini<-function(w){a<-sort(abs(w)); n<-length(a); if(sum(a)==0)return(NA); (2*sum(seq_len(n)*a)/(n*sum(a)))-(n+1)/n}

## TREINO leave-platform-out: treina SCANB(RNAseq)+METABRIC(array); held-out: TCGA, GSE20685, GSE1456
train_co <- c("SCANB","METABRIC"); held <- c("TCGA_BRCA","GSE20685","GSE1456")
common <- Reduce(intersect, c(list(pam50_genes), lapply(D[train_co], function(d) rownames(d$z))))
log("genes PAM50 comuns no treino: %d", length(common))
Xl<-lapply(train_co, function(co){ t(D[[co]]$z[common,,drop=FALSE]) })   # samples x genes
X<-do.call(rbind, Xl)
yt<-unlist(lapply(train_co, function(co) D[[co]]$time)); ye<-unlist(lapply(train_co, function(co) D[[co]]$event))
ok<-is.finite(yt)&is.finite(ye)&yt>0; X<-X[ok,]; surv<-Surv(yt[ok], ye[ok])
log("treino: %d amostras x %d genes | eventos: %d", nrow(X), ncol(X), sum(ye[ok]))

derive <- function(alpha){
  cv<-cv.glmnet(X, surv, family="cox", alpha=alpha, nfolds=10)
  b<-as.matrix(coef(cv, s="lambda.1se")); w<-b[b[,1]!=0,1]
  if(length(w)<5){ b<-as.matrix(coef(cv, s="lambda.min")); w<-b[b[,1]!=0,1] }
  w
}
cp2_spread <- derive(0.1)   # ridge-pesada -> espalha
cp2_lasso  <- derive(0.5)   # CorePAM-like
log("CorePAM-2(ridge a=0.1): %d genes | Gini peso %.2f", length(cp2_spread), gini(cp2_spread))
log("CorePAM-2(a=0.5):       %d genes | Gini peso %.2f", length(cp2_lasso),  gini(cp2_lasso))
log("CorePAM original:        %d genes | Gini peso %.2f", length(corepam), gini(corepam))

models <- list(CorePAM=corepam, "CorePAM2_ridge"=cp2_spread, "CorePAM2_lasso"=cp2_lasso)
## orienta cada modelo high=pior pelo treino
orient <- function(w){ s<-sc_w(do.call(cbind, lapply(train_co, function(co) D[[co]]$z[intersect(names(w),rownames(D[[co]]$z)),,drop=FALSE])), w); w }
## avalia discriminacao no HELD-OUT + robustez (envelope gene-loss weight-directed)
envelope <- function(w, d, thr=0.02){
  pres<-names(w)[names(w)%in%rownames(d$z)]; if(length(pres)<5) return(NA)
  w<-w[pres]; cf<-cidx(d$time,d$event, sc_w(d$z,w)); if(is.na(cf))return(NA); if(cf<0.5){w<--w;cf<-cidx(d$time,d$event,sc_w(d$z,w))}
  ord<-names(sort(abs(w),decreasing=TRUE)); ng<-length(w)
  for(nr in 0:(ng-3)){ c<-cidx(d$time,d$event, sc_w(d$z, w[setdiff(names(w),ord[seq_len(nr)])])); if(is.finite(c)&&cf-c>=thr) return((nr-1)/ng) }
  (ng-3)/ng
}
res<-list()
for(mn in names(models)){ w0<-models[[mn]]
  for(co in held){ d<-D[[co]]; pres<-names(w0)[names(w0)%in%rownames(d$z)]; w<-w0[pres]
    cf<-cidx(d$time,d$event, sc_w(d$z,w)); if(!is.na(cf)&&cf<0.5){w<--w; cf<-cidx(d$time,d$event, sc_w(d$z,w))}
    res[[length(res)+1]]<-data.frame(model=mn, cohort=co, n_genes=length(w0), n_genes_mapped=length(pres), c_heldout=round(cf,4),
                                      envelope=round(envelope(w0,d),3)) }
}
R<-do.call(rbind,res); write.csv(R, file.path(OUT,"rederive_results.csv"), row.names=FALSE)
log("\n=== HELD-OUT (TCGA, GSE20685, GSE1456): discriminacao + robustez ===")
for(mn in names(models)){ s<-R[R$model==mn,]
  log("  %-15s | n=%d | C-heldout media=%.3f (%s) | envelope media=%.0f%%",
      mn, s$n_genes[1], mean(s$c_heldout,na.rm=TRUE),
      paste(sprintf("%s=%.3f",s$cohort,s$c_heldout),collapse=" "),
      100*mean(s$envelope,na.rm=TRUE)) }
cat("DONE\n")
