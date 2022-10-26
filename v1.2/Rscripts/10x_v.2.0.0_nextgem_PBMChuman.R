########################################################################################
#
# INSTALL REQUIRED PACKAGES IF NEEDED
# 
# N.B. This may require the installation of local libraries. Please check the README
# file of the project for a list of required packages
#
# ######################################################################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.15")

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils',
                       'HDF5Array', 'terra', 'ggrastr'))



if (!requireNamespace("devtools", quietly = TRUE)) 
  install.packages("devtools", dependencies = c("Depends"))
devtools::install_github('cole-trapnell-lab/monocle3')

BiocManager::install(c("Gviz", "GenomicRanges", "rtracklayer"))
devtools::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")

if (!requireNamespace("tidyverse", quietly = TRUE)) 
  install.packages("tidyverse", dependencies = c("Depends"))

if (!requireNamespace("aricode", quietly = TRUE)) 
  install.packages("aricode", dependencies = c("Depends"))

if (!requireNamespace("Seurat", quietly = TRUE)) 
  install.packages("Seurat", dependencies = c("Depends"))

if (!requireNamespace("Seurat", quietly = TRUE)) 
  install.packages("R.utils", dependencies = c("Depends"))

if (!requireNamespace("SeuratWrappers", quietly = TRUE)) 
  remotes::install_github('satijalab/seurat-wrappers')

if (!requireNamespace("Signac", quietly = TRUE)) 
  install.packages("Signac", dependencies = c("Depends"))

if (!requireNamespace("readsparse", quietly = TRUE)) 
  install.packages("readsparse", dependencies = c("Depends"))

if (!requireNamespace("leiden", quietly = TRUE)) 
  install.packages("leiden", dependencies = c("Depends"))

if (!requireNamespace("hdf5r", quietly = TRUE)) 
  install.packages("hdf5r", dependencies = c("Depends"))

install.packages('BiocManager', dependencies = c("Depends"))
BiocManager::install()
BiocManager::install('multtest')
install.packages('Seurat')

if (!requireNamespace("Seurat", quietly = TRUE)) 
  install.packages("Seurat", dependencies = TRUE)

if (!requireNamespace("Seurat", quietly = TRUE)) 
  install.packages("R.utils", dependencies = c("Depends"))

if (!requireNamespace("SeuratWrappers", quietly = TRUE)) 
  remotes::install_github('satijalab/seurat-wrappers')

if (!requireNamespace("data.table", quietly = TRUE)) 
  install.packages("data.table", dependencies = c("Depends"))

if (!requireNamespace("cellranger", quietly = TRUE)) 
  install.packages("cellranger", dependencies = c("Depends"))

if (!requireNamespace("readxl", quietly = TRUE)) 
  install.packages("readxl", dependencies = c("Depends"))

if (!requireNamespace("sctransform", quietly = TRUE)) 
  install.packages("sctransform", dependencies = c("Depends"))

if (!requireNamespace("class", quietly = TRUE)) 
  install.packages("class", dependencies = c("Depends"))

if (!requireNamespace("GenomeInfoDb", quietly = TRUE)) 
  install.packages("GenomeInfoDb", dependencies = c("Depends"))


#######################################################################################
#
# LOAD REQUIRED PACKAGES
# 
#######################################################################################

library(cicero)
library(data.table)
library(Matrix)
library(proxy)
library(reshape2)
library(readsparse)
library(aricode)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Signac)
library(Seurat)
library(GenomeInfoDb)

#######################################################################################
#
# 10k V2.0.0 PBMC PROCESSING
# 
#######################################################################################
if (!(dir.exists("../TMPResults"))){
  dir.create("../TMPResults")
}

if (!(dir.exists("../TMPResults/classifications"))){
  dir.create("../TMPResults/classifications")
}
if (!(dir.exists("../TMPResults/classifications/10X_V2_PBMC"))){
  dir.create("../TMPResults/classifications/10X_V2_PBMC")
}

if (!(dir.exists("../TMPResults/GAM"))){
  dir.create("../TMPResults/GAM")
}
if (!(dir.exists("../TMPResults/GAM/10X_V2_PBMC"))){
  dir.create("../TMPResults/GAM/10X_V2_PBMC")
}

if (!(dir.exists("../TMPResults/metrics"))){
  dir.create("../TMPResults/metrics")
}
if (!(dir.exists("../TMPResults/metrics/10X_V2_PBMC"))){
  dir.create("../TMPResults/metrics/10X_V2_PBMC")
}

if (!(dir.exists("../TMPResults/Robjects"))){
  dir.create("../TMPResults/Robjects")
}
if (!(dir.exists("../TMPResults/Robjects/10X_V2_PBMC"))){
  dir.create("../TMPResults/Robjects/10X_V2_PBMC")
}

####### MATRIX LOADING AND ATAC ANALYSIS ########
v2_pbmc_peaks <- read.table("../TMPDATA/10X_V2_PBMC/atac_pbmc_5k_nextgem_peaks.bed")
v2_pbmc_peaks <- paste0(v2_pbmc_peaks$V1, "_", v2_pbmc_peaks$V2, "_", v2_pbmc_peaks$V3)
v2_pbmc_matrix <- Read10X_h5("../TMPDATA/10X_V2_PBMC/atac_pbmc_5k_nextgem_filtered_peak_bc_matrix.h5")
rownames(v2_pbmc_matrix) <- v2_pbmc_peaks

genome_ref = read.table("../DATA/Gene_2022/Genomes/hg38/hg38.p13.chrom.sizes.txt")
genome_ref <- genome_ref[1:24,]

hg38 <- genome_ref[1:24,]
hg38 <- Seqinfo(hg38$V1, seqlengths= hg38$V2)
hg38@genome[] <- "hg38"


v2_pbmc_matrix@x[v2_pbmc_matrix@x >0] <- 1

processed_ATAC_cds <- new_cell_data_set(v2_pbmc_matrix)
processed_ATAC_cds <- processed_ATAC_cds[,Matrix::colSums(exprs(processed_ATAC_cds)) != 0]
processed_ATAC_cds <- processed_ATAC_cds[Matrix::rowSums(exprs(processed_ATAC_cds)) != 0,]

processed_ATAC_cds <- detect_genes(processed_ATAC_cds)
processed_ATAC_cds <- estimate_size_factors(processed_ATAC_cds)
processed_ATAC_cds <- preprocess_cds(processed_ATAC_cds,method ="LSI")
processed_ATAC_cds <- reduce_dimension(processed_ATAC_cds, reduction_method = 'UMAP',  preprocess_method = "LSI")
processed_ATAC_cds <- cluster_cells(processed_ATAC_cds)
#plot_cells(processed_ATAC_cds)
saveRDS(processed_ATAC_cds, "../TMPResults/Robject/10X_V2_PBMC/processed_ATAC_cds")


####### CO-ACCESSIBILITY ########
genome_ref = read.table("../DATA/Gene_2022/Genomes/hg38/hg38.p13.chrom.sizes.txt")
genome_ref <- genome_ref[1:24,]
umap_coords <- reducedDims(processed_ATAC_cds)$UMAP
cicero_cds <- make_cicero_cds(processed_ATAC_cds, reduced_coordinates = umap_coords)
connection_table <- run_cicero(cicero_cds, genome_ref)

saveRDS(connection_table, "../TMPResults/Robject/10X_V2_PBMC/connection_table")
#connection_table <- readRDS("../TMPResults/connection_table_10x_v2.0.0")


con_val <- connection_table[connection_table$coaccess > 0,]
con_val <- con_val[!is.na(con_val$coaccess),]
coaccess <- signif(mean(con_val$coaccess), digits = 2)

####### SIGNAC ######

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)

PBMC_assay <- CreateChromatinAssay(
  counts = v1_pbmc_matrix,
  sep = c("_", "_"),
  genome = hg38,
  min.cells = 1
)

PBMC <- CreateSeuratObject(
  counts = PBMC_assay,
  assay = 'peaks',
  project = 'ATAC'
)

saveRDS(PBMC, "../TMPResults/Robject/10X_V2_PBMC/PBMC")

####### LABELING PEAKS #####

labeled_peaks <- read.csv("../TMPresults/labeled_peaks/10X_V2_PBMC/cCRE_labeled_peaks.tsv", sep = "\t")
#labeled_peaks_multi <- read.csv("C:/Users/loren/IWBBIO/data/labeled peaks/classifiedPeaks_multiCols.csv", sep = "\t")

nmax <- max(stringr::str_count(labeled_peaks$encodeCcreCombined_hg38_ucscLabel, "\t")) + 1

labeled_peaks <- separate(labeled_peaks, col = encodeCcreCombined_hg38_ucscLabel, sep = "\t", into = paste0("RegFunc", seq_len(nmax)))

labeled_peaks$site_names <- paste0(labeled_peaks$X.chrom, "_", labeled_peaks$chromStart, "_", labeled_peaks$chromEnd)

labeled_peaks <- labeled_peaks[labeled_peaks$site_names %in% rownames(fData(processed_ATAC_cds)),]
labeled_peaks <- labeled_peaks[!duplicated(labeled_peaks),]

saveRDS(labeled_peaks, "../TMPResults/Robject/10X_V2_PBMC/labeled_peaks")
saveRDS(nmax, "../TMPResults/Robject/10X_V2_PBMC/nmax")


