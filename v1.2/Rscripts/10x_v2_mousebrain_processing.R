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
# 10k PBMC Multiome Controller PROCESSING
# 
#######################################################################################
if (!(dir.exists("../TMPResults"))){
  dir.create("../TMPResults")
}

if (!(dir.exists("../TMPResults/classifications"))){
  dir.create("../TMPResults/classifications")
}
if (!(dir.exists("../TMPResults/classifications/10x_v2_mousebrain"))){
  dir.create("../TMPResults/classifications/10x_v2_mousebrain")
}

if (!(dir.exists("../TMPResults/GAM"))){
  dir.create("../TMPResults/GAM")
}
if (!(dir.exists("../TMPResults/GAM/10x_v2_mousebrain"))){
  dir.create("../TMPResults/GAM/10x_v2_mousebrain")
}

if (!(dir.exists("../TMPResults/metrics"))){
  dir.create("../TMPResults/metrics")
}
if (!(dir.exists("../TMPResults/metrics/10x_v2_mousebrain"))){
  dir.create("../TMPResults/metrics/10x_v2_mousebrain")
}

if (!(dir.exists("../TMPResults/Robjects"))){
  dir.create("../TMPResults/Robjects")
}
if (!(dir.exists("../TMPResults/Robjects/10x_v2_mousebrain"))){
  dir.create("../TMPResults/Robjects/10x_v2_mousebrain")
}

####### MATRIX LOADING ########
matrix <- readMM("../TMPDATA/10x_v2_mousebrain/matrix.mtx")
cells <- read.table("../TMPDATA/10x_v2_mousebrain/barcodes.tsv")
features <- read.delim("../TMPDATA/10x_v2_mousebrain/peaks.bed", header=FALSE)

peaks <- paste0(features$V1,"_",features$V2, "_", features$V3)

row.names(matrix) <- peaks
colnames(matrix) <- cells$V1

matrix@x[matrix@x > 0] <- 1
processed_ATAC_cds <- new_cell_data_set(matrix)


processed_ATAC_cds <- detect_genes(processed_ATAC_cds)
processed_ATAC_cds <- estimate_size_factors(processed_ATAC_cds)
processed_ATAC_cds <- preprocess_cds(processed_ATAC_cds,method ="LSI")
processed_ATAC_cds <- reduce_dimension(processed_ATAC_cds, reduction_method = 'UMAP',  preprocess_method = "LSI")
processed_ATAC_cds <- cluster_cells(processed_ATAC_cds)

#plot_cells(processed_ATAC_cds)
saveRDS(processed_ATAC_cds, "../TMPResults/Robjects/10x_v2_mousebrain/processed_ATAC_cds")

####### CO-ACCESSIBILITY ########
#data("mouse.mm9.genome")
mm10 <- read.table("../DATA/Gene_2022/Genomes/mm10/mm10.chrom.sizes.txt")

umap_coords <- reducedDims(processed_ATAC_cds)$UMAP
cicero_cds <- make_cicero_cds(processed_ATAC_cds, reduced_coordinates = umap_coords)


connection_table <- run_cicero(cicero_cds, mm10.1, sample_num = 100)
saveRDS(connection_table, "../TMPResults/Robjects/10x_v2_mousebrain/connection_table")
#connection_table <- readRDS("../TMPResults/connection_table_10x_v2_mousebrain")

con_val <- connection_table[connection_table$coaccess > 0,]
con_val <- con_val[!is.na(con_val$coaccess),]
coaccess <- signif(mean(con_val$coaccess), digits = 2)

mm10.1 <- mm10[1:21,]
mm10.1 <- Seqinfo(mm10.1$V1, seqlengths= mm10.1$V2)

mm10.1@genome[] <- "mm10"

####### SIGNAC ############

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)


brain_assay <- CreateChromatinAssay(
  counts = matrix,
  sep = c("_", "_"),
  genome = mm10.1,
  min.cells = 1
)

brain <- CreateSeuratObject(
  counts = brain_assay,
  assay = 'peaks',
  project = 'ATAC'
)

saveRDS(brain, "../TMPResults/Robjects/10x_v2_mousebrain/brain")

####### LABELING PEAKS########

labeled_peaks <- read.csv("../TMPresults/labeled_peaks/10x_v2_mousebrain/cCRE_labeled_peaks.tsv", sep = "\t")

nmax <- max(stringr::str_count(labeled_peaks$encodeCcreCombined_ucscLabel, "\t")) + 1

labeled_peaks <- separate(labeled_peaks, col = encodeCcreCombined_ucscLabel, sep = "\t", into = paste0("RegFunc", seq_len(nmax)))


labeled_peaks$site_names <- paste0(labeled_peaks$X.chrom, "_", labeled_peaks$chromStart, "_", labeled_peaks$chromEnd)

labeled_peaks <- labeled_peaks[labeled_peaks$site_names %in% rownames(fData(processed_ATAC_cds)),]
labeled_peaks <- labeled_peaks[!duplicated(labeled_peaks),]

saveRDS(labeled_peaks, "../TMPResults/Robjects/10x_v2_mousebrain/labeled_peaks")
saveRDS(nmax, "../TMPResults/Robjects/10x_v2_mousebrain/nmax")

####### CICERO COMPUTATION ###############

chr2acc <- read.csv("../DATA/Gene_2022/Genomes/mm10/refseq gene annotation/chr2acc.txt", sep = "\t")
gene_anno <- rtracklayer::readGFF("../DATA/Gene_2022/Genomes/mm10/refseq gene annotation/GCF_000001635.26_GRCm38.p6_genomic.gtf.gz")
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

processed_ATAC_cds <- annotate_cds_by_site(processed_ATAC_cds, gene_annotation_sub)

tail(fData(processed_ATAC_cds))




unnorm_ga <- build_gene_activity_matrix(processed_ATAC_cds, connection_table)
unnorm_ga <- unnorm_ga[!Matrix::rowSums(unnorm_ga) == 0, 
                       !Matrix::colSums(unnorm_ga) == 0]
num_genes <- pData(processed_ATAC_cds)$num_genes_expressed
names(num_genes) <- row.names(pData(processed_ATAC_cds))

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
cds_cicero = cluster_cells(cds_cicero, resolution=1e-3)

#plot_cells(cds_cicero)

class <- as.data.frame(processed_ATAC_cds@clusters@listData[["UMAP"]][["clusters"]])
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
write.table(c_gam, file='../TMPResults/GAM/10x_v2_mousebrain/cicero.tsv', quote=FALSE, sep='\t', col.names = NA)
write.table(class2, file='../TMPResults/classifications/10x_v2_mousebrain/cicero_classification.tsv', quote=FALSE, sep='\t', col.names = NA)
write.table(class, file='../TMPResults/classifications/10x_v2_mousebrain/atac_classification.tsv', quote=FALSE, sep='\t', col.names = NA)

####### GENE SCORING COMPUTATION ################


library(GenomicRanges)
library(SummarizedExperiment)
library(data.table)
library(dplyr)
library(BuenColors)
library(Matrix)


chr2acc <- read.csv("../DATA/Gene_2022/Genomes/mm10/refseq gene annotation/chr2acc.txt", sep = "\t")
gene_anno <- rtracklayer::readGFF("../DATA/Gene_2022/Genomes/mm10/refseq gene annotation/GCF_000001635.26_GRCm38.p6_genomic.gtf.gz")
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

adf <- data.frame(do.call(rbind,strsplit(rownames(matrix),'_')),stringsAsFactors = FALSE)
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
fm_genescoring <- (t %*% matrix)
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

class <- as.data.frame(processed_ATAC_cds@clusters@listData[["UMAP"]][["clusters"]])
colnames(class) <- "CLASS"

class3 <- as.data.frame(cds_gs@clusters@listData[["UMAP"]][["clusters"]])

colnames(class3) <- "CLASS"
cds_gs@colData@listData[["ATAC"]] <- class$CLASS

#plot_cells(cds_gs, color_cells_by = "ATAC", label_groups_by_cluster = FALSE)

ARI <- ARI(class$CLASS, class3$CLASS)
AMI<- AMI(class$CLASS, class3$CLASS)
ARI
AMI

tsv <- data.matrix(exprs(cds_gs))
write.table(tsv, file='../TMPResults/GAM/10x_v2_mousebrain/genescoring.tsv', quote=FALSE, sep='\t', col.names = NA)
write.table(class3, file='../TMPResults/classifications/10x_v2_mousebrain/genescoring_classification.tsv', quote=FALSE, sep='\t', col.names = NA)






