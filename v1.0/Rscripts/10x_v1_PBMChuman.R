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

### 10x_V1_PBMChuman is the folder for the dataset elements; it is the folder where to put the dtasets of human PBMC
#  peaks.bed -> peaks file
#  matrix.h5 - > scATAC-seq matrix
#  conversion_error.txt
#
# follows loading of the dataset

v1_pbmc_peaks <- read.table("../10x_V1_PBMChuman/peaks.bed")
v1_pbmc_peaks <- paste0(v1_pbmc_peaks$V1, "_", v1_pbmc_peaks$V2, "_", v1_pbmc_peaks$V3)
v1_pbmc_matrix <- Read10X_h5("../10x_V1_PBMChuman/matrix.h5")
rownames(v1_pbmc_matrix) <- v1_pbmc_peaks

##### only for conversion errors
conversion_error <- read.table("../conversion_error.txt")
conversion_error <- paste0(conversion_error$V1,"_",conversion_error$V2, "_", conversion_error$V3)

diff <- setdiff(v1_pbmc_peaks, conversion_error)

v1_pbmc_matrix <- v1_pbmc_matrix[diff,]

genome_ref = read.table("../DATA/IWBBIO_2022/Genomes/hg38/hg38.p13.chrom.sizes.txt")
genome_ref <- genome_ref[1:24,]


hg38peaks_v1_pbmc <- read.table("../DATA/10x_V1_PBMChuman/hglft_genome_327f4_5c430.bed")
new_peaks <- paste0(hg38peaks_v1_pbmc$V1,"_",hg38peaks_v1_pbmc$V2, "_", hg38peaks_v1_pbmc$V3)

new_new_peaks <- hg38peaks_v1_pbmc[hg38peaks_v1_pbmc$V1 %in% genome_ref$V1,]
new_new_peaks <- paste0(new_new_peaks$V1,"_",new_new_peaks$V2, "_", new_new_peaks$V3)
new_new_peaks <- new_new_peaks[!duplicated(new_new_peaks)]

diff2 <- setdiff(new_peaks, new_new_peaks)
diff3 <- setdiff(new_peaks, diff2)
rownames(v1_pbmc_matrix) <- new_peaks
v1_pbmc_matrix <- v1_pbmc_matrix[diff3,]
rownames(v1_pbmc_matrix) <- new_new_peaks
#####


v1_pbmc_matrix@x[v1_pbmc_matrix@x >0] <- 1

input_cds <- new_cell_data_set(v1_pbmc_matrix)
input_cds <- input_cds[,Matrix::colSums(exprs(input_cds)) != 0]
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,]

input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds,method ="LSI")
input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP',  preprocess_method = "LSI")
input_cds <- cluster_cells(input_cds)
#plot_cells(input_cds)


umap_coords <- reducedDims(input_cds)$UMAP
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)
conns <- run_cicero(cicero_cds, genome_ref)


saveRDS(conns, "../TMPResults/conns_synthetic")

con_val <- conns[conns$coaccess > 0,]
con_val <- con_val[!is.na(con_val$coaccess),]
coaccess <- signif(mean(con_val$coaccess), digits = 2)

###############




gene_anno <- rtracklayer::readGFF("../DATA/IWBBIO_2022/Genomes/hg38/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz")
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
gene_annotation_sub <- gene_annotation_sub[,c("chromosome", "start", "end", "symbol")]

# Rename the gene symbol column to "gene"
names(gene_annotation_sub)[4] <- "gene"

input_cds <- annotate_cds_by_site(input_cds, gene_annotation_sub)

tail(fData(input_cds))

unnorm_ga <- build_gene_activity_matrix(input_cds, conns, coaccess_cutoff = coaccess)
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
cds_cicero = cluster_cells(cds_cicero, resolution=6e-3)

#plot_cells(cds_cicero)

#class <- as.data.frame(input_cds@clusters@listData[["UMAP"]][["clusters"]])
class2 <- as.data.frame(cds_cicero@clusters@listData[["UMAP"]][["clusters"]])
#colnames(class) <- "CLASS"
colnames(class2) <- "CLASS"
cds_cicero@colData@listData[["label"]] <- synthetic_metadata$label

#plot_cells(cds_cicero, color_cells_by = "label", label_groups_by_cluster = FALSE)

#marker_test_res_rna <- top_markers(cds_cicero)
#plot_cells(cds_cicero, genes = "Spi1")

c_gam <- data.matrix(exprs(cds_cicero))
write.table(c_gam, file='../TMPResults/cicero_gam.tsv', quote=FALSE, sep='\t', col.names = NA)
write.table(class2, file='../TMPResults/classification/cicero_classification.tsv', quote=FALSE, sep='\t', col.names = NA)


ARI(synthetic_metadata$label, class2$CLASS)
AMI(synthetic_metadata$label, class2$CLASS)

####


library(GenomicRanges)
library(SummarizedExperiment)
library(data.table)
library(dplyr)
library(BuenColors)
library(Matrix)

gene_anno <- rtracklayer::readGFF("../DATA/IWBBIO_2022/Genomes/hg38/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz")
gene_anno$chromosome <- gene_anno$seqid
#gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
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

#ov <- findOverlaps(atacgranges, tss_idx)
#options(repr.plot.width=3, repr.plot.height=3)
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

label <- as.data.frame(input_cds@colData@listData)
class3 <- as.data.frame(cds_gs@clusters@listData[["UMAP"]][["clusters"]])
colnames(label)[1] <- "CLASS"
colnames(class3) <- "CLASS"
cds_gs@colData@listData[["ATAC"]] <- label$CLASS

#plot_cells(cds_gs, color_cells_by = "ATAC", label_groups_by_cluster = FALSE)

#plot_cells(cds_gs)

gs_gam <- data.matrix(exprs(cds_gs))
write.table(c_gam, file='../TMPResults/genescoring_gam.tsv', quote=FALSE, sep='\t', col.names = NA)
write.table(class3, file='../TMPResults/classification/genescoring_classification.tsv', quote=FALSE, sep='\t', col.names = NA)


ARI(label$CLASS, class3$CLASS)
AMI(label$CLASS, class3$CLASS)

#########

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)

PBMC_assay <- CreateChromatinAssay(
  counts = synthetic_matrix,
  sep = c("_", "_"),
  genome = hg38,
  min.cells = 1
)

PBMC <- CreateSeuratObject(
  counts = PBMC_assay,
  assay = 'peaks',
  project = 'ATAC'
)




########

labeled_peaks <- read.csv("../TMPResults/leabeled_peaks/encodeCcreCombined_hg38_ucscLabel_classifiedPeaks.csv", sep = "\t")
#labeled_peaks_multi <- read.csv("C:/Users/loren/IWBBIO/data/labeled peaks/classifiedPeaks_multiCols.csv", sep = "\t")

nmax <- max(stringr::str_count(labeled_peaks$encodeCcreCombined_hg38_ucscLabel, "\t")) + 1

labeled_peaks <- separate(labeled_peaks, col = encodeCcreCombined_hg38_ucscLabel, sep = "\t", into = paste0("RegFunc", seq_len(nmax)))

labeled_peaks$site_names <- paste0(labeled_peaks$X.chrom, "_", labeled_peaks$chromStart, "_", labeled_peaks$chromEnd)

labeled_peaks <- labeled_peaks[labeled_peaks$site_names %in% rownames(fData(input_cds)),]
labeled_peaks <- labeled_peaks[!duplicated(labeled_peaks),]

peaks_multi_info <- fData(input_cds)
#rownames(peaks_multi_info) <- gsub(":", "_",rownames(peaks_multi_info))
#rownames(peaks_multi_info) <- gsub("-", "_",rownames(peaks_multi_info))
peaks_multi_info$site_name <- rownames(peaks_multi_info)
peaks_multi_info$is.prom <- FALSE
peaks_multi_info$is.enhD <- FALSE
peaks_multi_info <- cbind(peaks_multi_info,labeled_peaks[7:486])



ppp <- peaks_multi_info
ppp <- as.data.frame(peaks_multi_info)
ppp[ppp == "prom"] <- TRUE

peaks_multi_info$is.prom <- apply(ppp, 1, any)
peaks_multi_info[is.na(peaks_multi_info$is.prom),]$is.prom <- FALSE



peaks_prom <- peaks_multi_info[peaks_multi_info$is.prom == "TRUE",]
#different_prom <- peaks_prom[is.na(peaks_prom$gene),]

prom_list <- rownames(peaks_prom)
non_prom_list <- rownames(peaks_multi_info)
non_prom_list <- setdiff(non_prom_list, prom_list)

non_prom_peaks <- peaks_multi_info[non_prom_list,]


#### refseq #####
#refseq_gene_anno <-  rtracklayer::readGFF("C:/Users/loren/OneDrive - Politecnico di Torino (1)/IWBBIO/DATA/genomes/hg19/hg19.ncbiRefSeq.gtf.gz")
refseq_gene_anno <-  rtracklayer::readGFF("../DATA/IWBBIO_2022/Genomes/hg38/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz")

chr2acc <- read.csv("../DATA/Genomes/hg38/chr2acc.txt", sep = "\t")


refseq_peaks_prom <- peaks_multi_info[peaks_multi_info$is.prom == "TRUE",]
refseq_prom_list <- rownames(refseq_peaks_prom)
#refseq_prom_list <- refseq_prom_list[!duplicated(refseq_prom_list)]

refseq_gene_anno <- refseq_gene_anno[refseq_gene_anno$seqid %in% chr2acc$Accession.version,]
refseq_gene_anno$seqid <- as.factor(as.character(refseq_gene_anno$seqid))
levels(refseq_gene_anno$seqid) <- chr2acc$X.Chromosome
refseq_gene_anno$seqid <- paste0("chr", refseq_gene_anno$seqid)
refseq_gene_anno <- refseq_gene_anno[refseq_gene_anno$type == "gene",]
refseq_gene_anno <- refseq_gene_anno[refseq_gene_anno$gene_biotype %in% c("protein_coding","lncRNA"),]

hg38 <- genome_ref[1:24,]
#hg38$V1 <- as.character.factor(hg38$V1)
hg38 <- Seqinfo(hg38$V1, seqlengths= hg38$V2)
hg38@genome[] <- "hg38"



refseq_GRanges <- makeGRangesFromDataFrame(refseq_gene_anno, seqinfo = hg38, seqnames.field = "seqid", keep.extra.columns = TRUE)

refseq_gene_closest_to_prom <- ClosestFeature(PBMC, gsub("_", "-", refseq_prom_list), annotation = refseq_GRanges)
refseq_gene_closest_to_prom <- refseq_gene_closest_to_prom[refseq_gene_closest_to_prom$distance <= 500,]

refseq_near_prom_list <- gsub("-","_",refseq_gene_closest_to_prom$query_region)


refseq_peaks_prom <-  refseq_peaks_prom[refseq_near_prom_list,]
refseq_peaks_prom$gene_anno <- refseq_gene_closest_to_prom$gene
refseq_peaks_prom$site_name <- rownames(refseq_peaks_prom)

refseq_promoter_peak_table <- refseq_peaks_prom[,c("gene_anno","site_name")]
#refseq_promoter_peak_table$site_name <- rownames(refseq_promoter_peak_table)
prom_gene_name <- levels(factor(refseq_promoter_peak_table$gene_anno))

refseq_promoter_gene_mat <-
  Matrix::sparseMatrix(j=as.numeric(factor(refseq_promoter_peak_table$site_name)),
                       i=as.numeric(factor(refseq_promoter_peak_table$gene_anno)),
                       x=1)

refseq_accessibility_mat <- exprs(input_cds)
#rownames(refseq_accessibility_mat) <- gsub(":", "_",rownames(refseq_accessibility_mat))
#rownames(refseq_accessibility_mat) <- gsub("-", "_",rownames(refseq_accessibility_mat))

refseq_promoter_activity_scores <- refseq_accessibility_mat[refseq_near_prom_list,, drop=FALSE]

colnames(refseq_promoter_gene_mat) = levels(factor(refseq_promoter_peak_table$site_name))
row.names(refseq_promoter_gene_mat) = levels(factor(refseq_promoter_peak_table$gene_anno))
refseq_promoter_gene_mat <- refseq_promoter_gene_mat[,row.names(refseq_promoter_activity_scores)]

refseq_first_gene_matrix2 <- refseq_promoter_gene_mat %*% refseq_promoter_activity_scores
refseq_first_gene_matrix <- refseq_promoter_gene_mat %*% refseq_promoter_activity_scores
refseq_first_gene_matrix@x[refseq_first_gene_matrix@x > 0] <- 1

######### intragenic #######

intragenetic_non_prom_peaks <- ClosestFeature(brain, gsub("_", "-", non_prom_list), refseq_GRanges)
intragenetic_non_prom_peaks <- intragenetic_non_prom_peaks[intragenetic_non_prom_peaks$distance == 0,]
intragenetic_non_prom_peaks_list <- gsub("-","_",intragenetic_non_prom_peaks$query_region)

intragenetic_peaks <-  non_prom_peaks[intragenetic_non_prom_peaks_list,]
intragenetic_peaks$gene_anno <- intragenetic_non_prom_peaks$gene
intragenetic_peaks$site_name <- rownames(intragenetic_peaks)

intragenetic_peak_table <- intragenetic_peaks[,c("gene_anno","site_name")]


peak_table <- rbind(refseq_promoter_peak_table, intragenetic_peak_table)
allgene_peaks_list <- rownames(peak_table)

refseq_allgene_gene_mat <-
  Matrix::sparseMatrix(j=as.numeric(factor(peak_table$site_name)),
                       i=as.numeric(factor(peak_table$gene_anno)),
                       x=1)
colnames(refseq_allgene_gene_mat) = levels(factor(peak_table$site_name))
row.names(refseq_allgene_gene_mat) = levels(factor(peak_table$gene_anno))
refseq_allgene_gene_mat <- refseq_allgene_gene_mat[prom_gene_name,]

refseq_allgene_accessibility_mat <- exprs(CDS_ATAC)
rownames(refseq_allgene_accessibility_mat) <- gsub(":", "_",rownames(refseq_allgene_accessibility_mat))
rownames(refseq_allgene_accessibility_mat) <- gsub("-", "_",rownames(refseq_allgene_accessibility_mat))

refseq_allgene_activity_scores <- refseq_allgene_accessibility_mat[allgene_peaks_list,, drop=FALSE]

refseq_allgene_gene_mat <- refseq_allgene_gene_mat[,row.names(refseq_allgene_activity_scores)]

refseq_allgene_gene_matrix <- refseq_allgene_gene_mat %*% refseq_allgene_activity_scores

#refseq_promoter_activity_scores <- refseq_accessibility_mat[refseq_near_prom_list,, drop=FALSE]


refseq_allgene_gene_matrix <- refseq_allgene_gene_matrix * refseq_first_gene_matrix

###################


accessibility_mat <- exprs(input_cds)

if ("dist" %in% colnames(conns) == FALSE) {
  Peak1_cols <- split_peak_names(conns$Peak1)
  Peak2_cols <- split_peak_names(conns$Peak2)
  Peak1_bp <- round((as.integer(Peak1_cols[,3]) +
                       as.integer(Peak1_cols[,2])) / 2)
  Peak2_bp <- round((as.integer(Peak2_cols[,3]) +
                       as.integer(Peak2_cols[,2])) / 2)
  conns$dist <- abs(Peak2_bp - Peak1_bp)
}

site_weights = NULL
if (is.null(site_weights)) {
  site_weights <- Matrix::rowMeans(accessibility_mat) /
    Matrix::rowMeans(accessibility_mat)
  site_weights[names(site_weights)] <- 1
}
site_names <- names(site_weights)
site_weights <- as(Matrix::Diagonal(x=as.numeric(site_weights)),
                   "sparseMatrix")
row.names(site_weights) <- site_names
colnames(site_weights) <- site_names




ppp <- peaks_multi_info
ppp <- as.data.frame(peaks_multi_info)
ppp[ppp == "enhD"] <- TRUE

peaks_multi_info$is.enhD <- FALSE
peaks_multi_info$is.enhD <- apply(ppp,1, any)
peaks_multi_info[is.na(peaks_multi_info$is.enhD),]$is.enhD <- FALSE

#ppp <- peaks_multi_info[peaks_multi_info$is.prom == "TRUE",]

enhD_peaks_list <- rownames(peaks_multi_info[peaks_multi_info$is.enhD == TRUE,])
enhD_peaks_list <- setdiff(enhD_peaks_list, prom_list)


promoter_peak_table <- refseq_peaks_prom[, c("site_name", "gene_anno")]


promoter_peak_table$site_name <- as.character(row.names(promoter_peak_table))
promoter_peak_table <- promoter_peak_table[!is.na(promoter_peak_table$gene_anno),]
promoter_peak_table <- promoter_peak_table[,c("site_name", "gene_anno")]
promoter_peak_table$gene_anno <- as.character(promoter_peak_table$gene_anno)

colnames(promoter_peak_table) <- c("peak", "gene")

dist_thresh=300000


prom_enhD <- conns[(conns$Peak1 %in%
                      promoter_peak_table$peak &
                      conns$Peak2 %in%
                      enhD_peaks_list) | (conns$Peak2 %in%
                                            promoter_peak_table$peak &
                                            conns$Peak1 %in%
                                            enhD_peaks_list),]
prom_enhD <- prom_enhD[prom_enhD$coaccess >= coaccess & prom_enhD$dist <= dist_thresh,]
prom_enhD <- prom_enhD[!duplicated(prom_enhD),]

prom_enhD <- prom_enhD[,c("Peak1", "Peak2", "coaccess")]
prom_enhD <- prom_enhD[!duplicated(prom_enhD),]

prom_enhD$Peak1 <- as.character(prom_enhD$Peak1)
prom_enhD$Peak2 <- as.character(prom_enhD$Peak2)


prom_enhD <- rbind(prom_enhD,
                   data.frame(Peak1=unique(promoter_peak_table$peak),
                              Peak2=unique(promoter_peak_table$peak),
                              coaccess=0))

prom_enhD_connectivity <- make_sparse_matrix(prom_enhD, x.name = "coaccess")
#distal_connectivity_matrix <- make_sparse_matrix(nonneg_cons, x.name="coaccess")




promoter_conn_matrix2 <-
  prom_enhD_connectivity[unique(promoter_peak_table$peak),]

# Get list of promoter and distal sites in accessibility mat
promoter_safe_sites2 <- intersect(rownames(promoter_conn_matrix2),
                                  row.names(accessibility_mat))
distal_safe_sites2 <- intersect(colnames(promoter_conn_matrix2),
                                row.names(accessibility_mat))
distal_safe_sites2 <- setdiff(distal_safe_sites2, promoter_safe_sites2)

# Get accessibility info for promoters
promoter_access_mat_in_cicero_map2 <- accessibility_mat[promoter_safe_sites2,, drop=FALSE]

# Get accessibility for distal sites
distal_activity_scores2 <- accessibility_mat[distal_safe_sites2,, drop=FALSE]

# Scale connectivity matrix by site_weights
scaled_site_weights2 <- site_weights[distal_safe_sites2,distal_safe_sites2, drop=FALSE]
total_linked_site_weights2 <- promoter_conn_matrix2[,distal_safe_sites2, drop=FALSE] %*%
  scaled_site_weights2
total_linked_site_weights2 <- 1/Matrix::rowSums(total_linked_site_weights2,
                                                na.rm=TRUE)
total_linked_site_weights2[is.finite(total_linked_site_weights2) == FALSE] <- 0
total_linked_site_weights2[is.na(total_linked_site_weights2)] <- 0
total_linked_site_weights2[is.nan(total_linked_site_weights2)] <- 0
total_linked_site_weights2 <- Matrix::Diagonal(x=total_linked_site_weights2)
scaled_site_weights2 <- total_linked_site_weights2 %*%
  promoter_conn_matrix2[,distal_safe_sites2, drop=FALSE] %*%
  scaled_site_weights2
scaled_site_weights2@x[scaled_site_weights2@x > 1] <- 1

# Multiply distal accessibility by site weights
distal_activity_scores2 <- scaled_site_weights2 %*% distal_activity_scores2

distal_activity_scores2 <-
  distal_activity_scores2[row.names(promoter_access_mat_in_cicero_map2),, drop=FALSE]

promoter_activity_scores2 <- refseq_promoter_activity_scores + distal_activity_scores2

promoter_gene_mat2 <-
  Matrix::sparseMatrix(j=as.numeric(factor(promoter_peak_table$peak)),
                       i=as.numeric(factor(promoter_peak_table$gene)),
                       x=1)
colnames(promoter_gene_mat2) = levels(factor(promoter_peak_table$peak))
row.names(promoter_gene_mat2) = levels(factor(promoter_peak_table$gene))
promoter_gene_mat2 <- promoter_gene_mat2[,row.names(promoter_activity_scores2)]
gene_activity_scores2 <- promoter_gene_mat2 %*% promoter_activity_scores2


##################


unnorm_first <- gene_activity_scores2[!Matrix::rowSums(gene_activity_scores2) == 0, 
                                      !Matrix::colSums(gene_activity_scores2) == 0]
num_genes <- pData(input_cds)$num_genes_expressed
names(num_genes) <- row.names(pData(input_cds))

#norm_first_gene_matrix <- normalize_gene_activities(unnorm_first, num_genes)

final_GAM <- norm_first_gene_matrix #+ refseq_first_gene_matrix2

###########
final_GAM <- unnorm_first



first_cell <- colnames(final_GAM)
#cicero_cell <- read.table(cicero_cell)
lenght1 <- length(first_cell)

first_gene <- row.names(final_GAM)
lenght2 <- length(first_gene)
c_c <- matrix(first_cell, nrow = lenght1, dimnames = list(first_cell,c("Cells")))
c_g<- matrix(first_gene, nrow = lenght2, dimnames = list(first_gene,c("gene_short_name")))


final_GAM_cds <-  suppressWarnings(new_cell_data_set(final_GAM, cell_metadata = c_c, gene_metadata = c_g))

final_GAM_cds <- detect_genes(final_GAM_cds)
final_GAM_cds <- estimate_size_factors(final_GAM_cds)
final_GAM_cds <- preprocess_cds(final_GAM_cds, method = "PCA")

final_GAM_cds <- reduce_dimension(final_GAM_cds, reduction_method = 'UMAP', 
                                  preprocess_method = "PCA")
final_GAM_cds = cluster_cells(final_GAM_cds, resolution= 2e-3)

#plot_cells(final_GAM_cds, cell_size = 0.8, group_label_size =5)


#class <- as.data.frame(CDS_ATAC@clusters@listData[["UMAP"]][["clusters"]])
#colnames(class) <- "CLASS"

final_GAM_cds@colData@listData[["label"]] <- synthetic_metadata$label

#classification$CellType
#final_GAM_cds@colData@listData[["ATAC"]] <- class$CLASS
#cell_type_cells <- row.names(subset(pData(final_GAM_cds), label !="Ambiguous" & label != "Unknown"))

#plot_cells(final_GAM_cds, color_cells_by = "ATAC", label_groups_by_cluster = FALSE ,cell_size = 1.2, group_label_size = 5, label_cell_groups = FALSE)


#plot_cells(final_GAM_cds, color_cells_by = "label", label_groups_by_cluster = FALSE , group_label_size = 5, label_cell_groups = FALSE)


final_GAM_first <- as.data.frame(final_GAM_cds@clusters@listData[["UMAP"]][["clusters"]])
colnames(final_GAM_first) <- "CLASS"

tsv <- data.matrix(exprs(final_GAM_cds))

write.table(tsv, file='../TMPResults/gagam.tsv', quote=FALSE, sep='\t', col.names = NA)
write.table(final_GAM_first, file='../TMPResults/gagam_classification.tsv', quote=FALSE, sep='\t', col.names = NA)


ARI(synthetic_metadata$label, final_GAM_first$CLASS)
AMI(synthetic_metadata$label, final_GAM_first$CLASS)



#########functions ####
split_peak_names <- function(inp) {
  out <- stringr::str_split_fixed(stringi::stri_reverse(inp), 
                                  ":|-|_", 3)
  out[,1] <- stringi::stri_reverse(out[,1])
  out[,2] <- stringi::stri_reverse(out[,2])
  out[,3] <- stringi::stri_reverse(out[,3])
  out[,c(3,2,1), drop=FALSE]
}
make_sparse_matrix <- function(data,
                               i.name = "Peak1",
                               j.name = "Peak2",
                               x.name = "value") {
  if(!i.name %in% names(data) |
     !j.name %in% names(data) |
     !x.name %in% names(data)) {
    stop('i.name, j.name, and x.name must be columns in data')
  }
  
  data$i <- as.character(data[,i.name])
  data$j <- as.character(data[,j.name])
  data$x <- data[,x.name]
  
  if(!class(data$x) %in%  c("numeric", "integer"))
    stop('x.name column must be numeric')
  
  peaks <- data.frame(Peak = unique(c(data$i, data$j)),
                      index = seq_len(length(unique(c(data$i, data$j)))))
  
  data <- data[,c("i", "j", "x")]
  
  data <- rbind(data, data.frame(i=peaks$Peak, j = peaks$Peak, x = 0))
  data <- data[!duplicated(data[,c("i", "j", "x")]),]
  data <- data.table::as.data.table(data)
  peaks <- data.table::as.data.table(peaks)
  data.table::setkey(data, "i")
  data.table::setkey(peaks, "Peak")
  data <- data[peaks]
  data.table::setkey(data, "j")
  data <- data[peaks]
  data <- as.data.frame(data)
  
  data <- data[,c("index", "i.index", "x")]
  data2 <- data
  names(data2) <- c("i.index", "index", "x")
  
  data <- rbind(data, data2)
  
  data <- data[!duplicated(data[,c("index", "i.index")]),]
  data <- data[data$index >= data$i.index,]
  
  sp_mat <- Matrix::sparseMatrix(i=as.numeric(data$index),
                                 j=as.numeric(data$i.index),
                                 x=data$x,
                                 symmetric = TRUE)
  
  colnames(sp_mat) <- peaks[order(peaks$index),]$Peak
  row.names(sp_mat) <- peaks[order(peaks$index),]$Peak
  return(sp_mat)
}




