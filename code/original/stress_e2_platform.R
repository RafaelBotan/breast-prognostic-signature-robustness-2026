# Part of the reproducibility repository for:
# "Robustness and transportability of breast-cancer prognostic gene signatures
#  under structured gene loss and platform change."
# Public cohorts only (GEO/TCGA/cBioPortal); no patient-level data are redistributed.
# Set COHORT_DIR to the folder with the processed cohort matrices (see README).

# E2 — transportabilidade de PLATAFORMA (deployment): intra-coorte vs frozen single-sample
# Trava principio 6: referencia congelada vem de coortes OUTRAS (leave-cohort-out), nunca da alvo.
# intra  = z-score dentro da coorte-alvo (modo pesquisa, melhor caso)
# frozen_within = referencia = media/DP por gene das OUTRAS coortes da MESMA plataforma (nested)
# frozen_cross  = referencia = SCANB (RNA-seq) aplicada ao alvo (choque de plataforma explicito)
suppressMessages({library(arrow); library(survival)})
base <- Sys.getenv("COHORT_DIR", "data/cohorts")
OUT  <- "results"
log <- function(...) cat(sprintf(...), "\n")
suppressMessages({library(genefu)})

## assinaturas (mesma construcao do E1)
corepam <- c(EXO1=0.2174,NAT1=-0.2052,BLVRA=-0.1848,ACTR3B=-0.1405,MIA=-0.1193,MYBL2=-0.1183,
  PTTG1=0.1140,MDM2=-0.1074,SFRP1=-0.0907,GPR160=-0.0841,FOXC1=0.0717,PHGDH=0.0620,MYC=-0.0424,
  ESR1=0.0396,KRT5=-0.0277,CXXC5=0.0266,PGR=-0.0262,KRT17=-0.0163,CENPF=0.0150,FGFR4=0.0110,
  BCL2=-0.0105,MLPH=-0.0099,ERBB2=-0.0082,GRB7=-0.0053)
data(sig.oncotypedx); data(sig.ggi); data(sig.gene70)
od<-sig.oncotypedx; od<-od[od$group!="reference",]; oncotype<-setNames(od$weight*ifelse(od$group=="estrogen",-1,1), od$symbol)
gg<-sig.ggi[!is.na(sig.ggi$HUGO.gene.symbol),]; ggi<-setNames(ifelse(gg$grade==max(gg$grade,na.rm=TRUE),1,-1), gg$HUGO.gene.symbol); ggi<-ggi[!duplicated(names(ggi))]
g70<-sig.gene70[!is.na(sig.gene70$HUGO.gene.symbol),]; mammaprint<-setNames(-g70$correlation, g70$HUGO.gene.symbol); mammaprint<-mammaprint[!duplicated(names(mammaprint))]
sigs <- list(CorePAM=corepam, OncotypeDX=oncotype, GGI=ggi, MammaPrint=mammaprint)

cohorts <- c("SCANB","TCGA_BRCA","METABRIC","GSE20685","GSE1456")
fam <- c(SCANB="rnaseq", TCGA_BRCA="rnaseq", METABRIC="array", GSE20685="array", GSE1456="array")

## carrega tudo (expressao crua gene x amostra + tempo/evento alinhados)
D <- list()
for(co in cohorts){
  ex<-as.data.frame(read_parquet(sprintf("%s/%s/expression_genelevel_preZ.parquet",base,co)))
  g<-ex$gene; ex$gene<-NULL; ex<-as.matrix(ex); rownames(ex)<-g
  ex<-ex[, !grepl("REPL$",colnames(ex)),drop=FALSE]
  cl<-as.data.frame(read_parquet(sprintf("%s/%s/clinical_FINAL.parquet",base,co)))
  key<-if(sum(colnames(ex)%in%cl$patient_id)>sum(colnames(ex)%in%cl$sample_id))"patient_id" else "sample_id"
  m<-match(colnames(ex), cl[[key]])
  D[[co]]<-list(ex=ex, time=cl$os_time_months[m], event=cl$os_event[m])
  log("carregada %s: %dx%d", co, nrow(ex), ncol(ex))
}
cidx<-function(time,event,sc){ok<-is.finite(time)&is.finite(event)&is.finite(sc); if(sum(ok)<20||length(unique(sc[ok]))<3)return(NA_real_); tryCatch(concordance(Surv(time[ok],event[ok])~sc[ok])$concordance,error=function(e)NA_real_)}
sc_w<-function(z,w){p<-names(w)[names(w)%in%rownames(z)]; if(length(p)<3)return(rep(NA_real_,ncol(z))); as.numeric(colSums(z[p,,drop=FALSE]*w[p])/sum(abs(w[p])))}

## referencia congelada (mean/sd por gene) de um conjunto de coortes
ref_params<-function(cohs){
  exl<-lapply(cohs, function(c) D[[c]]$ex)
  common<-Reduce(intersect, lapply(exl, rownames))
  M<-do.call(cbind, lapply(exl, function(e) e[common,,drop=FALSE]))
  list(mean=rowMeans(M,na.rm=TRUE), sd=apply(M,1,sd,na.rm=TRUE))
}
z_frozen<-function(ex, rp){ g<-intersect(rownames(ex), names(rp$mean)); z<-(ex[g,,drop=FALSE]-rp$mean[g])/rp$sd[g]; z[is.finite(rowSums(z)),,drop=FALSE] }

rows<-list(); add<-function(...){rows[[length(rows)+1]]<<-data.frame(..., stringsAsFactors=FALSE)}
for(sg in names(sigs)){ w<-sigs[[sg]]
  for(co in cohorts){ d<-D[[co]]
    z_intra<-t(scale(t(d$ex))); z_intra<-z_intra[is.finite(rowSums(z_intra)),,drop=FALSE]
    c_in<-cidx(d$time,d$event, sc_w(z_intra,w)); flip<- !is.na(c_in)&&c_in<0.5
    ww<- if(flip) -w else w; if(flip) c_in<-cidx(d$time,d$event, sc_w(z_intra,ww))
    add(sig=sg, target=co, plat=fam[[co]], mode="intra", cindex=round(c_in,4), n=ncol(d$ex))
    # frozen within (leave-cohort-out, mesma plataforma)
    others<-setdiff(names(fam)[fam==fam[[co]]], co)
    if(length(others)>=1){ rp<-ref_params(others); zf<-z_frozen(d$ex,rp)
      add(sig=sg, target=co, plat=fam[[co]], mode="frozen_within", cindex=round(cidx(d$time,d$event, sc_w(zf,ww)),4), n=ncol(d$ex)) }
    # frozen cross (referencia SCANB RNA-seq)
    if(co!="SCANB"){ rp<-ref_params("SCANB"); zf<-z_frozen(d$ex,rp)
      add(sig=sg, target=co, plat=fam[[co]], mode="frozen_cross", cindex=round(cidx(d$time,d$event, sc_w(zf,ww)),4), n=ncol(d$ex)) }
  }
}
res<-do.call(rbind, rows)
write.csv(res, file.path(OUT,"e2_transport.csv"), row.names=FALSE)

## resumo: gap intra - frozen, por modo, por assinatura
log("\n=== C-index por modo (media; principio 6: ver distribuicao no CSV) ===")
for(sg in names(sigs)){
  s<-res[res$sig==sg,]
  ci<-mean(s$cindex[s$mode=="intra"],na.rm=TRUE)
  fw<-mean(s$cindex[s$mode=="frozen_within"],na.rm=TRUE)
  fc<-mean(s$cindex[s$mode=="frozen_cross"],na.rm=TRUE)
  log("  %-11s intra=%.3f | frozen_within=%.3f (gap %.3f) | frozen_cross=%.3f (gap %.3f)",
      sg, ci, fw, ci-fw, fc, ci-fc)
}
log("\n=== detalhe CorePAM por coorte ===")
s<-res[res$sig=="CorePAM",]
for(co in cohorts){ r<-s[s$target==co,]
  g<-function(md){v<-r$cindex[r$mode==md]; if(length(v))sprintf("%.3f",v) else "-"}
  log("  %-10s [%s] intra=%s within=%s cross=%s", co, fam[[co]], g("intra"), g("frozen_within"), g("frozen_cross")) }
cat("DONE\n")
