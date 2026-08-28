suppressPackageStartupMessages({
  library(arrow)
  library(GEOquery)
  library(genefu)
  library(survival)
})

cache_dir <- Sys.getenv(
  "REVISION_CACHE_DIR",
  file.path(tempdir(), "cci_gse_cache")
)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
gse <- getGEO(
  "GSE7390",
  GSEMatrix = TRUE,
  AnnotGPL = TRUE,
  destdir = cache_dir
)[[1]]
p <- pData(gse)
print(grep("dmfs|distant|metast", colnames(p), ignore.case = TRUE, value = TRUE))
characteristic_cols <- grep("characteristics", colnames(p), ignore.case = TRUE)
for (j in characteristic_cols) {
  values <- as.character(p[[j]])
  if (any(grepl("dmfs|distant|metast", values, ignore.case = TRUE))) {
    cat("\nCOLUMN", colnames(p)[j], "\n")
    print(head(values, 12))
    print(summary(as.numeric(sub(".*: *", "", values))))
  }
}
