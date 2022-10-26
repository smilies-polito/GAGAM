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
# 10k PBMC Multiome Chromium X PROCESSING
# 
#######################################################################################
if (!(dir.exists("../TMPResults"))){
  dir.create("../TMPResults")
}

if (!(dir.exists("../TMPResults/classifications"))){
  dir.create("../TMPResults/classifications")
}
if (!(dir.exists("../TMPResults/classifications/10x_PBMC_Multiome_ChromiumX"))){
  dir.create("../TMPResults/classifications/10x_PBMC_Multiome_ChromiumX")
}

if (!(dir.exists("../TMPResults/GAM"))){
  dir.create("../TMPResults/GAM")
}
if (!(dir.exists("../TMPResults/GAM/10x_PBMC_Multiome_ChromiumX"))){
  dir.create("../TMPResults/GAM/10x_PBMC_Multiome_ChromiumX")
}

if (!(dir.exists("../TMPResults/metrics"))){
  dir.create("../TMPResults/metrics")
}
if (!(dir.exists("../TMPResults/metrics/10x_PBMC_Multiome_ChromiumX"))){
  dir.create("../TMPResults/metrics/10x_PBMC_Multiome_ChromiumX")
}

if (!(dir.exists("../TMPResults/Robjects"))){
  dir.create("../TMPResults/Robjects")
}
if (!(dir.exists("../TMPResults/Robjects/10x_PBMC_Multiome_ChromiumX"))){
  dir.create("../TMPResults/Robjects/10x_PBMC_Multiome_ChromiumX")
}
####### MATRIX LOADING AND CONVERSION ERROR CORRECTION ########

matrix <- readMM("../DATA/10x_PBMC_Multiome_ChromiumX/matrix.mtx")
cells <- read.table("../DATA/10x_PBMC_Multiome_ChromiumX/barcodes.tsv")
features <- read.delim("../DATA/10x_PBMC_Multiome_ChromiumX/features.tsv", header=FALSE)

#division of the fetures between genes and peaks
genes <- features[features$V3 == "Gene Expression",]
colnames(genes)[2] <- "gene_short_name"
peaks <- features[features$V3 == "Peaks",]
peaks$V1 <- gsub("-","_",peaks$V1)
peaks$V1 <- gsub(":","_",peaks$V1)
peaks$V2 <- gsub("-","_",peaks$V2)
peaks$V2 <- gsub(":","_",peaks$V2)
#peaks <- peaks[,4:6]
#write.table(peaks, file='10x_multiome_X_peaks.csv', quote=FALSE, sep='\t', col.names = NA)

row.names(matrix) <- features$V2
colnames(matrix) <- cells$V1

#creation of the two matrices, accordingly to the features
RNA_matrix <- matrix[genes$gene_short_name,]
ATAC_matrix <- matrix[peaks$V2,]
ATAC_matrix@Dimnames[[1]] <- gsub("-","_",ATAC_matrix@Dimnames[[1]])
ATAC_matrix@Dimnames[[1]] <- gsub(":","_",ATAC_matrix@Dimnames[[1]])
####### RNA ANALYSIS #######
CDS_RNA <- new_cell_data_set(RNA_matrix)
#adding gene names as a metadata for rows, useful for later functions
rowData(CDS_RNA)$gene_short_name <- genes$gene_short_name
CDS_RNA <- detect_genes(CDS_RNA)
CDS_RNA <- estimate_size_factors(CDS_RNA)
#preprocessing consisting in normalization, scaling, and dimansional reduction (both LSI and PCA)
CDS_RNA <- preprocess_cds(CDS_RNA, method = "LSI")
CDS_RNA <- preprocess_cds(CDS_RNA, method = "PCA")
#non-linear dimensional reduction UMAP, to visualize the cells in 2D
CDS_RNA <- reduce_dimension(CDS_RNA, reduction_method = 'UMAP', 
                            preprocess_method = "LSI")
#clustering of the cells, based on the UMAP
CDS_RNA <- cluster_cells(CDS_RNA, resolution=0.6e-3)
#plot_cells(CDS_RNA,reduction_method = 'UMAP', group_label_size = 5) + ggtitle("RNA CLUSTERING")
saveRDS(CDS_RNA,"../TMPResults/Robjects/10x_PBMC_Multiome_ChromiumX/CDS_RNA")

marker_test_res_rna <- top_markers(CDS_RNA)
top_specific_markers_rna <- marker_test_res_rna %>%
  filter(fraction_expressing >= 0.10) %>%
  filter(specificity >= 0.15) %>%
  group_by(cell_group) %>%
  top_n(1, pseudo_R2)
top_specific_marker_ids <- unique(top_specific_markers_rna %>% pull(gene_id))
#plot_genes_by_group(CDS_RNA,
#                    top_specific_marker_ids,
 #                   group_cells_by="cluster",
  #                  ordering_type="maximal_on_diag")

####### ATAC ANALYSIS #######
##binarization of the matrix
ATAC_matrix@x[ATAC_matrix@x > 0] <- 1
#Creation of the CDS object for the ATAC data
CDS_ATAC <- new_cell_data_set(ATAC_matrix)
rowData(CDS_ATAC)$gene_short_name <- peaks$V2
#the process and function are totally analogous to  before
CDS_ATAC <- detect_genes(CDS_ATAC)
CDS_ATAC <- estimate_size_factors(CDS_ATAC)
CDS_ATAC <- preprocess_cds(CDS_ATAC, method = "LSI")
#CDS_ATAC <- preprocess_cds(CDS_ATAC, method = "PCA")
CDS_ATAC <- reduce_dimension(CDS_ATAC, reduction_method = 'UMAP', 
                             preprocess_method = "LSI")
CDS_ATAC <- cluster_cells(CDS_ATAC, resolution=0.8e-3)
#plot the results based on ATAC data alone
#plot_cells(CDS_ATAC,reduction_method = 'UMAP', group_label_size = 5)#+ ggtitle("ATAC CLUSTERING")

input_cds <- CDS_ATAC
#rm(CDS_ATAC)
saveRDS(input_cds, "../TMPResults/Robjects/10x_PBMC_Multiome_ChromiumX/input_cds")

####### CO-ACCESSIBILITY ########
genome_ref = read.table("../DATA/Gene_2022/Genomes/hg38/hg38.p13.chrom.sizes.txt")
genome_ref <- genome_ref[1:24,]

umap_coords <- reducedDims(input_cds)$UMAP
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)
conns <- run_cicero(cicero_cds, genome_ref)

saveRDS(conns, "../TMPResults/Robjects/10x_PBMC_Multiome_ChromiumX/conns")
#conns <- readRDS("../TMPResults/conns_10k_Multiome_X")

con_val <- conns[conns$coaccess > 0,]
con_val <- con_val[!is.na(con_val$coaccess),]
coaccess <- signif(mean(con_val$coaccess), digits = 2)

####### LABELING PEAKS ###########
labeled_peaks <- read.csv("../TMPresults/labeled_peaks/10x_PBMC_Multiome_ChromiumX/cCRE_labeled_peaks.tsv", sep = "\t")
nmax <- max(stringr::str_count(labeled_peaks$encodeCcreCombined_hg38_ucscLabel, "\t")) + 1

labeled_peaks <- separate(labeled_peaks, col = encodeCcreCombined_hg38_ucscLabel, sep = "\t", into = paste0("RegFunc", seq_len(nmax)))

labeled_peaks$site_names <- paste0(labeled_peaks$X.chrom, "_", labeled_peaks$chromStart, "_", labeled_peaks$chromEnd)

labeled_peaks <- labeled_peaks[labeled_peaks$site_names %in% rownames(fData(input_cds)),]
labeled_peaks <- labeled_peaks[!duplicated(labeled_peaks),]

saveRDS(labeled_peaks, "../TMPResults/Robjects/10x_PBMC_Multiome_ChromiumX/labeled_peaks")
saveRDS(nmax, "../TMPResults/Robjects/10x_PBMC_Multiome_ChromiumX/nmax")

####### SIGNAC #######


PBMC_assay <- CreateChromatinAssay(
  counts = ATAC_matrix,
  sep = c("_", "_"),
  genome = hg38,
  min.cells = 1
)

PBMC <- CreateSeuratObject(
  counts = PBMC_assay,
  assay = 'peaks',
  project = 'ATAC'
)

saveRDS(PBMC, "../TMPResults/Robjects/10x_PBMC_Multiome_ChromiumX/PBMC")

####### CICERO COMPUTATION #######
chr2acc <- read.csv("../DATA/Gene_2022/Genomes/hg38/chr2acc.txt", sep = "\t")
gene_anno <- rtracklayer::readGFF("../DATA/Gene_2022/Genomes/hg38/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz")
gene_anno <- gene_anno[gene_anno$seqid %in% chr2acc$Accession.version,]
gene_anno$seqid <- as.factor(as.character(gene_anno$seqid))
levels(gene_anno$seqid) <- chr2acc$X.Chromosome
gene_anno$seqid <- paste0("chr", gene_anno$seqid)


gene_anno$chromosome <- gene_anno$seqid
#gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene



pos <- subset(gene_anno, strand == "+")
pos <- pos[order(pos$start),] 
# remove all but the first exons per transcript
pos <- pos[!duplicated(pos$transcript),] 
# make a 1 base pair marker of the TSS
pos$end <- pos$start + 1 

neg <- subset(gene_anno, strand == "-")
neg <- neg[order(neg$start, decreasing = TRUE),] 
# remove all but the first exons per transcript
neg <- neg[!duplicated(neg$transcript),] 
neg$start <- neg$end - 1

gene_annotation_sub <- rbind(pos, neg)


# Make a subset of the TSS annotation columns containing just the coordinates 
# and the gene name
gene_annotation_sub <- gene_annotation_sub[,c("seqid", "start", "end", "gene_id")]


# Rename the gene symbol column to "gene"
names(gene_annotation_sub)[4] <- "gene"

input_cds <- annotate_cds_by_site(input_cds, gene_annotation_sub)

tail(fData(input_cds))




unnorm_ga <- build_gene_activity_matrix(input_cds, conns)
unnorm_ga <- unnorm_ga[!Matrix::rowSums(unnorm_ga) == 0, 
                       !Matrix::colSums(unnorm_ga) == 0]
num_genes <- pData(input_cds)$num_genes_expressed
names(num_genes) <- row.names(pData(input_cds))

cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)

cicero_cell <- colnames(cicero_gene_activities)
#cicero_cell <- read.table(cicero_cell)
lenght1 <- length(cicero_cell)

cicero_gene <- row.names(cicero_gene_activities)
lenght2 <- length(cicero_gene)
c_c <- matrix(cicero_cell, nrow = lenght1, dimnames = list(cicero_cell,c("Cells")))
c_g<- matrix(cicero_gene, nrow = lenght2, dimnames = list(cicero_gene,c("gene_short_name")))


## processing GAM with Cicero
cds_cicero <-  suppressWarnings(new_cell_data_set(cicero_gene_activities, cell_metadata = c_c, gene_metadata = c_g))

cds_cicero <- detect_genes(cds_cicero)
cds_cicero <- estimate_size_factors(cds_cicero)
cds_cicero <- preprocess_cds(cds_cicero, method = "PCA")

cds_cicero <- reduce_dimension(cds_cicero, reduction_method = 'UMAP', 
                               preprocess_method = "PCA")
cds_cicero = cluster_cells(cds_cicero, resolution=0.8e-3)

#plot_cells(cds_cicero)

class <- as.data.frame(input_cds@clusters@listData[["UMAP"]][["clusters"]])
class2 <- as.data.frame(cds_cicero@clusters@listData[["UMAP"]][["clusters"]])
colnames(class) <- "CLASS"
colnames(class2) <- "CLASS"
cds_cicero@colData@listData[["CLASS"]] <- class$CLASS

#plot_cells(cds_cicero, color_cells_by = "CLASS", label_groups_by_cluster = FALSE)

#marker_test_res_rna <- top_markers(cds_cicero)
#plot_cells(cds_cicero, genes = "Spi1")

ARI <- ARI(class$CLASS, class2$CLASS)
AMI <- AMI(class$CLASS, class2$CLASS)
ARI
AMI

#saveRDS(cds_cicero, "cicero_gam_multiome")

c_gam <- data.matrix(exprs(cds_cicero))
write.table(c_gam, file='../TMPResults/GAM/10x_PBMC_Multiome_ChromiumX/cicero.tsv', quote=FALSE, sep='\t', col.names = NA)
write.table(class2, file='../TMPResults/classifications/10x_PBMC_Multiome_ChromiumX/cicero_classification.tsv', quote=FALSE, sep='\t', col.names = NA)
write.table(class, file='../TMPResults/classifications/10x_PBMC_Multiome_ChromiumX/atac_classification.tsv', quote=FALSE, sep='\t', col.names = NA)

####### GENE SCORING COMPUTATION ################


library(GenomicRanges)
library(SummarizedExperiment)
library(data.table)
library(dplyr)
library(BuenColors)
library(Matrix)


chr2acc <- read.csv("../DATA/Gene_2022/Genomes/hg38/chr2acc.txt", sep = "\t")
gene_anno <- rtracklayer::readGFF("../DATA/Gene_2022/Genomes/hg38/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz")
gene_anno <- gene_anno[gene_anno$seqid %in% chr2acc$Accession.version,]
gene_anno$seqid <- as.factor(as.character(gene_anno$seqid))
levels(gene_anno$seqid) <- chr2acc$X.Chromosome
gene_anno$seqid <- paste0("chr", gene_anno$seqid)

gene_anno$chromosome <- gene_anno$seqid
#gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene


gdf <- gene_anno[,c("chromosome", "start", "end", "symbol", "strand", "type")]
gdf <- gdf[gdf$type == "gene",]
#gdf <- read.table("C:/Users/loren/OneDrive - Politecnico di Torino (1)/IWBBIO/DATA/genomes/hg19/hg19-tss.bed.txt", stringsAsFactors = FALSE)

dim(gdf)
gdf[1:3,1:3]
colnames(gdf) <- c("V1", "V2", "V3", "V4", "V5", "V6")
tss <- data.frame(chr = gdf$V1, gene = gdf$V4, stringsAsFactors = FALSE)
tss$tss <-  ifelse(gdf$V5 == "+", gdf$V3, gdf$V2)
tss$start <- ifelse(tss$tss - 50000 > 0, tss$tss - 50000, 0)
tss$stop <- tss$tss + 50000
tss_idx <- makeGRangesFromDataFrame(tss, keep.extra.columns = TRUE)

adf <- data.frame(do.call(rbind,strsplit(rownames(ATAC_matrix),'_')),stringsAsFactors = FALSE)
colnames(adf) <- c("chr", "start", "end")
adf$start <- as.integer(adf$start)
adf$end <- as.integer(adf$end)
dim(adf)

adf$mp <- (adf$start + adf$end)/2
atacgranges <- makeGRangesFromDataFrame(adf, start.field = "mp", end.field = "mp")

ov <- findOverlaps(atacgranges, tss_idx)
options(repr.plot.width=3, repr.plot.height=3)
# plot a histogram showing peaks per gene
#qplot(table(subjectHits(ov)), binwidth = 1) + theme(plot.subtitle = element_text(vjust = 1), 
#                                                    plot.caption = element_text(vjust = 1)) +
#  labs(title = "Histogram of peaks per gene",  x = "Peaks / gene", y="Frequency") + pretty_plot()

dist <- abs(mcols(tss_idx)$tss[subjectHits(ov)] - start(atacgranges)[queryHits(ov)])
exp_dist_model <- exp(-1*dist/5000)


m <- Matrix::sparseMatrix(i = c(queryHits(ov), length(atacgranges)),
                          j = c(subjectHits(ov), length(tss_idx)),
                          x = c(exp_dist_model,0))

colnames(m) <- gdf$V4 # gene name
m <- m[,which(Matrix::colSums(m) != 0)]

#counts <- data.matrix(ATAC_matrix)
t <- t(m)
fm_genescoring <- (t %*% ATAC_matrix)
fm_genescoring <- fm_genescoring[!duplicated(rownames(fm_genescoring)),]

fm_genescoring[1:3, 1:3]
gene_s <- as.data.frame(rownames(fm_genescoring))
colnames(gene_s)[1] <- "gene_short_name"
rownames(gene_s) <- gene_s$gene_short_name

cds_gs <- new_cell_data_set(fm_genescoring, gene_metadata = gene_s)
cds_gs <- detect_genes(cds_gs)
cds_gs <- estimate_size_factors(cds_gs)
cds_gs <- preprocess_cds(cds_gs, method = "PCA")

cds_gs <- reduce_dimension(cds_gs, reduction_method = 'UMAP', 
                           preprocess_method = "PCA")
cds_gs = cluster_cells(cds_gs, resolution=0.8e-3)

#plot_cells(cds_gs)

class <- as.data.frame(input_cds@clusters@listData[["UMAP"]][["clusters"]])
colnames(class) <- "CLASS"

class3 <- as.data.frame(cds_gs@clusters@listData[["UMAP"]][["clusters"]])

colnames(class3) <- "CLASS"
cds_gs@colData@listData[["ATAC"]] <- class$CLASS

#plot_cells(cds_gs, color_cells_by = "ATAC", label_groups_by_cluster = FALSE)

ARI <- ARI(class$CLASS, class3$CLASS)
AMI <- AMI(class$CLASS, class3$CLASS)
ARI
AMI

tsv <- data.matrix(exprs(cds_gs))
write.table(tsv, file='../TMPResults/GAM/10x_PBMC_Multiome_ChromiumX/genescoring.tsv', quote=FALSE, sep='\t', col.names = NA)
write.table(class3, file='../TMPResults/classifications/10x_PBMC_Multiome_ChromiumX/genescoring_classification.tsv', quote=FALSE, sep='\t', col.names = NA)



