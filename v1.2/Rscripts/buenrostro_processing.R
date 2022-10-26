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
# Buenrostro2018 PROCESSING
# 
#######################################################################################
if (!(dir.exists("../TMPResults"))){
  dir.create("../TMPResults")
}

if (!(dir.exists("../TMPResults/classifications"))){
  dir.create("../TMPResults/classifications")
}
if (!(dir.exists("../TMPResults/classifications/buenrostro"))){
  dir.create("../TMPResults/classifications/buenrostro")
}

if (!(dir.exists("../TMPResults/GAM"))){
  dir.create("../TMPResults/GAM")
}
if (!(dir.exists("../TMPResults/GAM/buenrostro"))){
  dir.create("../TMPResults/GAM/buenrostro")
}

if (!(dir.exists("../TMPResults/metrics"))){
  dir.create("../TMPResults/metrics")
}
if (!(dir.exists("../TMPResults/metrics/buenrostro"))){
  dir.create("../TMPResults/metrics/buenrostro")
}

if (!(dir.exists("../TMPResults/Robjects"))){
  dir.create("../TMPResults/Robjects")
}
if (!(dir.exists("../TMPResults/Robjects/buenrostro"))){
  dir.create("../TMPResults/Robjects/buenrostro")
}

####### MATRIX LOADING AND CONVERSION ERROR CORRECTION ########
conversion_error <- read.table("../DATA/Gene_2022/Genomes/conversion_error.txt")
conversion_error <- paste0(conversion_error$V1,"_",conversion_error$V2, "_", conversion_error$V3)

orig_peaks <- read.table("../DATA/buenrostro/GSE96769_PeakFile_20160207.csv.bed")
orig_peaks <- paste0(orig_peaks$V1,"_",orig_peaks$V2, "_", orig_peaks$V3)

diff <- setdiff(orig_peaks, conversion_error)

buenrostro_matrix <- read.table("../DATA/buenrostro/GSE96769_scATACseq_counts.txt")
buenrostro_cells<- read.csv2("../DATA/buenrostro/cells.txt",sep = ";", header = FALSE,)
cells <- t(buenrostro_cells)
cells <- as.data.frame(cells)
cells <- cells$V1
mat <- Matrix::sparseMatrix(i=as.numeric((buenrostro_matrix$V1)),
                     j=as.numeric((buenrostro_matrix$V2)),
                     x=buenrostro_matrix$V3)
colnames(mat) <- cells
rownames(mat) <- orig_peaks
mat <- mat[diff,]
hg38peaks_buenrostro <- read.table("../DATA/buenrostro/hglft_genome_241cf_9a53c0.bed")
new_peaks <- paste0(hg38peaks_buenrostro$V1,"_",hg38peaks_buenrostro$V2, "_", hg38peaks_buenrostro$V3)
new_new_peaks <- hg38peaks_buenrostro[hg38peaks_buenrostro$V1 %in% genome_ref$V1,]
new_new_peaks <- paste0(new_new_peaks$V1,"_",new_new_peaks$V2, "_", new_new_peaks$V3)
diff2 <- setdiff(new_peaks, new_new_peaks)
diff3 <- setdiff(new_peaks, diff2)
rownames(mat) <- new_peaks
mat <- mat[diff3,]
rownames(mat) <- new_new_peaks
metadata <- read.table('../DATA/buenrostro/metadata.tsv.txt',
                       header = TRUE,
                       stringsAsFactors=FALSE,quote="",row.names=1)

mat <- mat[,rownames(metadata)]



####### CICERO COMPUTATION #######

processed_ATAC_cds <- new_cell_data_set(mat, cell_metadata = metadata)
processed_ATAC_cds <- processed_ATAC_cds[,Matrix::colSums(exprs(processed_ATAC_cds)) != 0]
processed_ATAC_cds <- processed_ATAC_cds[Matrix::rowSums(exprs(processed_ATAC_cds)) != 0,]

processed_ATAC_cds <- detect_genes(processed_ATAC_cds)
processed_ATAC_cds <- estimate_size_factors(processed_ATAC_cds)
processed_ATAC_cds <- preprocess_cds(processed_ATAC_cds,method ="LSI")
processed_ATAC_cds <- reduce_dimension(processed_ATAC_cds, reduction_method = 'UMAP',  preprocess_method = "LSI")
processed_ATAC_cds <- cluster_cells(processed_ATAC_cds, resolution = 3e-3)
#plot_cells(processed_ATAC_cds)
#plot_cells(processed_ATAC_cds, color_cells_by = "label", label_groups_by_cluster = FALSE, label_cell_groups = FALSE)

processed_ATAC_cds_first <- as.data.frame(processed_ATAC_cds@clusters@listData[["UMAP"]][["clusters"]])
colnames(processed_ATAC_cds_first) <- "CLASS"

#ARI(metadata$label, processed_ATAC_cds_first$CLASS)
#AMI(metadata$label, processed_ATAC_cds_first$CLASS)

saveRDS(processed_ATAC_cds, "../TMPResults/Robjects/buenrostro/processed_ATAC_cds")

####### CO-ACCESSIBILITY ########
#data("human.hg19.genome")
genome_ref = read.table("../DATA/Gene_2022/Genomes/hg38/hg38.p13.chrom.sizes.txt")
genome_ref <- genome_ref[1:24,]


umap_coords <- reducedDims(processed_ATAC_cds)$UMAP
cicero_cds <- make_cicero_cds(processed_ATAC_cds, reduced_coordinates = umap_coords)
connection_table <- run_cicero(cicero_cds, genome_ref)

saveRDS(connection_table, "../TMPResults/Robjects/buenrostro/connection_table")
#connection_table <- readRDS("../TMPResults/connection_table_buenrostro")


con_val <- connection_table[connection_table$coaccess > 0,]
con_val <- con_val[!is.na(con_val$coaccess),]
coaccess <- signif(mean(con_val$coaccess), digits = 2)

chr2acc <- read.csv("../DATA/Gene_2022/Genomes/hg38/chr2acc.txt", sep = "\t")
gene_anno <- rtracklayer::readGFF("../DATA/Gene_2022/Genomes/hg38/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz")
gene_anno <- gene_anno[gene_anno$seqid %in% chr2acc$Accession.version,]
gene_anno$seqid <- as.factor(as.character(gene_anno$seqid))
levels(gene_anno$seqid) <- chr2acc$X.Chromosome
gene_anno$seqid <- paste0("chr", gene_anno$seqid)


gene_anno$chromosome <- gene_anno$seqid
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

unnorm_ga <- build_gene_activity_matrix(processed_ATAC_cds, connection_table, coaccess_cutoff = coaccess)
unnorm_ga <- unnorm_ga[!Matrix::rowSums(unnorm_ga) == 0, 
                       !Matrix::colSums(unnorm_ga) == 0]
num_genes <- pData(processed_ATAC_cds)$num_genes_expressed
names(num_genes) <- row.names(pData(processed_ATAC_cds))

cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)

cicero_cell <- colnames(cicero_gene_activities)
lenght1 <- length(cicero_cell)

cicero_gene <- row.names(cicero_gene_activities)
lenght2 <- length(cicero_gene)
c_c <- matrix(cicero_cell, nrow = lenght1, dimnames = list(cicero_cell,c("Cells")))
c_g<- matrix(cicero_gene, nrow = lenght2, dimnames = list(cicero_gene,c("gene_short_name")))

####### CICERO PROCESSING #######
cds_cicero <-  suppressWarnings(new_cell_data_set(cicero_gene_activities, cell_metadata = c_c, gene_metadata = c_g))

cds_cicero <- detect_genes(cds_cicero)
cds_cicero <- estimate_size_factors(cds_cicero)
cds_cicero <- preprocess_cds(cds_cicero, method = "PCA")

cds_cicero <- reduce_dimension(cds_cicero, reduction_method = 'UMAP', 
                               preprocess_method = "PCA")
cds_cicero = cluster_cells(cds_cicero, resolution=6e-3)

#plot_cells(cds_cicero)

class2 <- as.data.frame(cds_cicero@clusters@listData[["UMAP"]][["clusters"]])
colnames(class2) <- "CLASS"
cds_cicero@colData@listData[["label"]] <- metadata$label

#plot_cells(cds_cicero, color_cells_by = "label", label_groups_by_cluster = FALSE)

ARI <- ARI(metadata$label, class2$CLASS)
AMI <- AMI(metadata$label, class2$CLASS)
ARI
AMI


#SAVE CICERO GAM AND CLUSTERING
c_gam <- data.matrix(exprs(cds_cicero))
write.table(c_gam, file='../TMPResults/GAM/buenrostro/cicero.tsv', quote=FALSE, sep='\t', col.names = NA)
write.table(class2, file='../TMPResults/classifications/buenrostro/cicero_classification.tsv', quote=FALSE, sep='\t', col.names = NA)
write.table(class, file='../TMPResults/classifications/buenrostro/atac_classification.tsv', quote=FALSE, sep='\t', col.names = NA)



####### GENE SCORING COMPUTATION #######

library(GenomicRanges)
library(SummarizedExperiment)
library(data.table)
library(dplyr)
library(BuenColors)
library(Matrix)

gene_anno <- rtracklayer::readGFF("../DATA/Gene_2022/Genomes/hg38/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz")
gene_anno$chromosome <- gene_anno$seqid
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene

gdf <- gene_anno[,c("chromosome", "start", "end", "symbol", "strand", "type")]

dim(gdf)
gdf[1:3,1:3]
colnames(gdf) <- c("V1", "V2", "V3", "V4", "V5", "V6")
tss <- data.frame(chr = gdf$V1, gene = gdf$V4, stringsAsFactors = FALSE)
tss$tss <-  ifelse(gdf$V5 == "+", gdf$V3, gdf$V2)
tss$start <- ifelse(tss$tss - 50000 > 0, tss$tss - 50000, 0)
tss$stop <- tss$tss + 50000
tss_idx <- makeGRangesFromDataFrame(tss, keep.extra.columns = TRUE)

adf <- data.frame(do.call(rbind,strsplit(rownames(mat),'_')),stringsAsFactors = FALSE)
colnames(adf) <- c("chr", "start", "end")
adf$start <- as.integer(adf$start)
adf$end <- as.integer(adf$end)
dim(adf)

adf$mp <- (adf$start + adf$end)/2
atacgranges <- makeGRangesFromDataFrame(adf, start.field = "mp", end.field = "mp")

ov <- findOverlaps(atacgranges, tss_idx)
options(repr.plot.width=3, repr.plot.height=3)

dist <- abs(mcols(tss_idx)$tss[subjectHits(ov)] - start(atacgranges)[queryHits(ov)])
exp_dist_model <- exp(-1*dist/5000)


m <- Matrix::sparseMatrix(i = c(queryHits(ov), length(atacgranges)),
                          j = c(subjectHits(ov), length(tss_idx)),
                          x = c(exp_dist_model,0))

colnames(m) <- gdf$V4 # gene name
m <- m[,which(Matrix::colSums(m) != 0)]

counts <- data.matrix(mat)
fm_genescoring <- data.matrix(t(m) %*% mat)
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
cds_gs = cluster_cells(cds_gs, resolution=3e-3)

label <- as.data.frame(processed_ATAC_cds@colData@listData)
class3 <- as.data.frame(cds_gs@clusters@listData[["UMAP"]][["clusters"]])
colnames(label)[1] <- "CLASS"
colnames(class3) <- "CLASS"
cds_gs@colData@listData[["ATAC"]] <- label$CLASS

#plot_cells(cds_gs, color_cells_by = "ATAC", label_groups_by_cluster = FALSE)

#plot_cells(cds_gs)

ARI <- ARI(label$CLASS, class3$CLASS)
AMI <- AMI(label$CLASS, class3$CLASS)
ARI
AMI
####### sIGNAC ########

exprs(processed_ATAC_cds)
PBMC_assay <- CreateChromatinAssay(
  counts = mat,
  sep = c("_", "_"),
  genome = hg38,
  min.cells = 1
)

PBMC <- CreateSeuratObject(
  counts = PBMC_assay,
  assay = 'peaks',
  project = 'ATAC'
)

saveRDS(PBMC, "../TMPResults/Robjects/buenrostro/brain")


####### LABELING PEAKS ###########

labeled_peaks <- read.csv("../TMPresults/labeled_peaks/buenrostro/cCRE_labeled_peaks.tsv", sep = "\t")

nmax <- max(stringr::str_count(labeled_peaks$encodeCcreCombined_hg38_ucscLabel, "\t")) + 1

labeled_peaks <- separate(labeled_peaks, col = encodeCcreCombined_hg38_ucscLabel, sep = "\t", into = paste0("RegFunc", seq_len(nmax)))

labeled_peaks$site_names <- paste0(labeled_peaks$X.chrom, "_", labeled_peaks$chromStart, "_", labeled_peaks$chromEnd)

labeled_peaks <- labeled_peaks[labeled_peaks$site_names %in% rownames(fData(processed_ATAC_cds)),]

saveRDS(labeled_peaks, "../TMPResults/Robjects/buenrostro/labeled_peaks")
saveRDS(nmax, "../TMPResults/Robjects/buenrostro/nmax")


#######################################################################################
#
# END
# 
#######################################################################################
