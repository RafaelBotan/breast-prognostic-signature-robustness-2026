# Part of the reproducibility repository for:
# "Robustness and transportability of breast-cancer prognostic gene signatures
#  under structured gene loss and platform change."
# Public cohorts only (GEO/TCGA/cBioPortal); no patient-level data are redistributed.
# Set COHORT_DIR to the folder with the processed cohort matrices (see README).

# Coorte CONFIRMATORIA (held-out NOVA, fora das 5) — GSE7390 (Desmedt/TRANSBIG)
# Roda a analise TRAVADA 1x: C-index das assinaturas + envelope (E1 weight-directed).
suppressMessages({library(GEOquery); library(survival)})
OUT <- "results"
log <- function(...) cat(sprintf(...), "\n")
options(timeout=600)

log("baixando GSE7390 (GEOquery)...")
gse <- tryCatch(getGEO("GSE7390", GSEMatrix=TRUE, AnnotGPL=TRUE, destdir=tempdir())[[1]],
                error=function(e){ log("ERRO download: %s", conditionMessage(e)); NULL })
if(is.null(gse)) quit(save="no", status=1)
ex <- exprs(gse); fd <- fData(gse); pd <- pData(gse)
log("expr: %d probes x %d amostras", nrow(ex), ncol(ex))
log("fData cols: %s", paste(head(colnames(fd),12),collapse=" | "))

## --- survival: achar tempo+evento (dmfs) ---
sym_col <- colnames(fd)[grep("symbol", colnames(fd), ignore.case=TRUE)][1]
log("coluna de simbolo: %s", sym_col)
# pData: procurar t.dmfs / e.dmfs em colunas ou em characteristics
findcol <- function(pat) { c<-colnames(pd)[grep(pat, colnames(pd), ignore.case=TRUE)]; if(length(c)) c[1] else NA }
tcol <- findcol("t.dmfs|t_dmfs|dmfs.time|time.dmfs"); ecol <- findcol("e.dmfs|e_dmfs|dmfs.event|event.dmfs")
parse_num <- function(x) as.numeric(gsub(".*: *","", as.character(x)))
if(is.na(tcol) || is.na(ecol)){
  # tentar characteristics_ch1.*
  ch <- pd[, grep("characteristics", colnames(pd)), drop=FALSE]
  flat <- apply(ch, 2, as.character)
  tcol2 <- which(apply(flat,2,function(col) any(grepl("t.dmfs|tdm", col, ignore.case=TRUE))))
  ecol2 <- which(apply(flat,2,function(col) any(grepl("e.dmfs|edm", col, ignore.case=TRUE))))
  if(length(tcol2)&&length(ecol2)){ time<-parse_num(ch[,tcol2[1]]); event<-parse_num(ch[,ecol2[1]]) }
  else { log("colunas de sobrevida (dmfs) NAO encontradas; cols pData: %s", paste(colnames(pd),collapse=" | ")); quit(save="no",status=2) }
} else { time<-parse_num(pd[[tcol]]); event<-parse_num(pd[[ecol]]) }
log("sobrevida: tempo[%s] evento[%s] | n eventos=%d | mediana tempo=%.0f", tcol, ecol, sum(event==1,na.rm=TRUE), median(time,na.rm=TRUE))

## --- probes -> simbolo (colapsa por maior variancia) ---
sym <- fd[[sym_col]]; sym <- sub(" ///.*","", sym)        # primeiro simbolo se multiplos
keep <- sym!="" & !is.na(sym)
ex<-ex[keep,]; sym<-sym[keep]
v <- apply(ex,1,var,na.rm=TRUE)
ord <- order(v, decreasing=TRUE)
ex<-ex[ord,]; sym<-sym[ord]
dup<-duplicated(sym); exg<-ex[!dup,]; rownames(exg)<-sym[!dup]
log("matriz gene-level: %d genes", nrow(exg))
z <- t(scale(t(exg))); z<-z[is.finite(rowSums(z)),,drop=FALSE]

## --- assinaturas (harmonizadas, = locked_stress_test) ---
suppressMessages(library(genefu)); data(sig.oncotypedx); data(sig.ggi); data(sig.gene70)
corepam <- c(EXO1=0.2174,NAT1=-0.2052,BLVRA=-0.1848,ACTR3B=-0.1405,MIA=-0.1193,MYBL2=-0.1183,PTTG1=0.1140,
  MDM2=-0.1074,SFRP1=-0.0907,GPR160=-0.0841,FOXC1=0.0717,PHGDH=0.0620,MYC=-0.0424,ESR1=0.0396,KRT5=-0.0277,
  CXXC5=0.0266,PGR=-0.0262,KRT17=-0.0163,CENPF=0.0150,FGFR4=0.0110,BCL2=-0.0105,MLPH=-0.0099,ERBB2=-0.0082,GRB7=-0.0053)
od<-sig.oncotypedx; od<-od[od$group!="reference",]; oncotype<-setNames(od$weight*ifelse(od$group=="estrogen",-1,1),od$symbol)
gg<-sig.ggi[!is.na(sig.ggi$HUGO.gene.symbol),]; ggi<-setNames(ifelse(gg$grade==max(gg$grade,na.rm=TRUE),1,-1),gg$HUGO.gene.symbol); ggi<-ggi[!duplicated(names(ggi))]
g70<-sig.gene70[!is.na(sig.gene70$HUGO.gene.symbol),]; mammaprint<-setNames(-g70$correlation,g70$HUGO.gene.symbol); mammaprint<-mammaprint[!duplicated(names(mammaprint))]
sigs<-list(CorePAM=corepam, OncotypeDX=oncotype, GGI=ggi, MammaPrint=mammaprint)
sc_w<-function(z,w){p<-names(w)[names(w)%in%rownames(z)]; if(length(p)<3)return(rep(NA_real_,ncol(z))); as.numeric(colSums(z[p,,drop=FALSE]*w[p])/sum(abs(w[p])))}
cidx<-function(t,e,s){ok<-is.finite(t)&is.finite(e)&is.finite(s); if(sum(ok)<20||length(unique(s[ok]))<3)return(NA); c<-tryCatch(concordance(Surv(t[ok],e[ok])~s[ok])$concordance,error=function(x)NA); c}

log("\n=== CONFIRMACAO GSE7390: C-index + envelope (weight-directed) ===")
rows<-list()
for(sg in names(sigs)){ w0<-sigs[[sg]]; pres<-names(w0)[names(w0)%in%rownames(z)]
  if(length(pres)<5){log("  %s: so %d genes -> pula",sg,length(pres)); next}
  w<-w0[pres]; cf<-cidx(time,event,sc_w(z,w)); if(!is.na(cf)&&cf<0.5){w<- -w; cf<-cidx(time,event,sc_w(z,w))}
  ng<-length(w); ord<-names(sort(abs(w),decreasing=TRUE)); env<-NA
  for(nr in 0:(ng-3)){ c<-cidx(time,event,sc_w(z,w[setdiff(names(w),ord[seq_len(nr)])])); if(is.finite(c)&&cf-c>=0.02){env<-(nr-1)/ng;break} }
  if(is.na(env)) env<-(ng-3)/ng
  log("  %-11s genes=%d/%d | C=%.3f | envelope=%.0f%%", sg, length(pres), length(w0), cf, 100*env)
  rows[[length(rows)+1]]<-data.frame(cohort="GSE7390",sig=sg,n_genes_present=length(pres),cindex=round(cf,4),envelope=round(env,3))
}
write.csv(do.call(rbind,rows), file.path(OUT,"confirmatory_GSE7390.csv"), row.names=FALSE)
cat("DONE\n")
