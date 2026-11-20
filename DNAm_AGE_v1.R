if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("methylclock", version = "devel")
BiocManager::install("IlluminaHumanMethylationEPICmanifest")
BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
BiocManager::install("methylclockData")
BiocManager::install("isglobal-brge/methylclock@devel")


library(methylclock)
library(sesame)
library(sesameData)
library(readxl)


excel_path <- "annotation.xlsx"
anno_df <- as.data.frame(read_excel(excel_path, sheet = 1))

data_dir <- "/file/path/idats"
basenames <- list.files(data_dir, pattern = "_Red\\.idat$", recursive = TRUE, full.names = TRUE)
basenames <- sub("_Red\\.idat$", "", basenames)


beta_list <- lapply(basenames, function(bn) {
  sdf <- readIDATpair(bn, manifest = anno_df)
  openSesame(sdf)
})
names(beta_list) <- basename(basenames)
betas_raw <- do.call(cbind, beta_list)

core_probes <- sub("_.*", "", rownames(betas_raw))


betas_clean <- rowsum(betas_raw, group = core_probes, na.rm = TRUE) / 
  rowsum((!is.na(betas_raw)) * 1, group = core_probes)


betas_clean <- as.matrix(betas_clean)
mode(betas_clean) <- "numeric"


clockCoverage(betas_clean)


res <- methylclock(
  betas_clean, 
  clocks = c("Horvath", "Levine", "TL"), 
  impute = "knn"
)

head(as.data.frame(res))

