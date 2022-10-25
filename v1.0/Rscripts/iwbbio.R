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



####loading datasets
# matrix.h5 -> scATAC-seq matrix
# peaks.bed -> peaks file


#binarized 
prova_10x <-Read10X_h5("../matrix.h5")
prova_10x@x[prova_10x@x > 0] <- 1
peak.info <- read.table("../peaks.bed")
peak.info$site_name <- paste(peak.info$V1, peak.info$V2, peak.info$V3, sep = "_")
row.names(peak.info) <- peak.info$site_name
row.names(prova_10x) <- row.names(peak.info)

input_cds <- new_cell_data_set(prova_10x, gene_metadata = peak.info)


input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds,method ="LSI")
input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP',  preprocess_method = "LSI")
input_cds <- cluster_cells(input_cds)

#plot_cells(input_cds)

mm10 <- read.table("../DATA/IWBBIO_2022/Genomes/mm10.chrom.sizes.txt")

umap_coords <- reducedDims(input_cds)$UMAP
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)


gene_ann_reg <- rtracklayer::readGFF("../DATA/IWBBIO_2022/Genomes/mus_musculus.GRCm38.Regulatory_Build.regulatory_features.20180516.gff.gz")
gene_ann_reg$seqid <- paste0("chr", gene_ann_reg$seqid)

gene_anno <-  rtracklayer::readGFF("../DATA/IWBBIO_2022/Genomes/gencode.vM17.annotation.gtf.gz")
gene_anno$chromosome <- gene_anno$seqid
#gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name

#data("human.hg19.genome")
genome_ref = human.hg19.genome

conns <- run_cicero(cicero_cds, mm10, sample_num = 100) 
saveRDS(conns, "../TMPResults/conns1")
conns <- readRDS("../TMPResults/conns1")

con_val <- conns[conns$coaccess > 0,]
con_val <- con_val[!is.na(con_val$coaccess),]
coaccess <- signif(mean(con_val$coaccess), digits = 2)


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
cds_cicero = cluster_cells(cds_cicero, resolution=0.95e-3)

#plot_cells(cds_cicero)

class <- as.data.frame(input_cds@clusters@listData[["UMAP"]][["clusters"]])
class2 <- as.data.frame(cds_cicero@clusters@listData[["UMAP"]][["clusters"]])
colnames(class) <- "CLASS"
colnames(class2) <- "CLASS"
cds_cicero@colData@listData[["ATAC"]] <- class$CLASS

#plot_cells(cds_cicero, color_cells_by = "ATAC")

marker_test_res_rna <- top_markers(cds_cicero)
#plot_cells(cds_cicero, genes = "Spi1")

ARI(class$CLASS, class2$CLASS)
AMI(class$CLASS, class2$CLASS)


###############gene scoring###################
library(GenomicRanges)
library(SummarizedExperiment)
library(data.table)
library(dplyr)
library(BuenColors)
library(Matrix)

gdf <- gene_anno[,c("chromosome", "start", "end", "symbol", "strand", "type")]
gdf <- gdf[gdf$type == "gene",]


colnames(gdf) <- c("V1", "V2", "V3", "V4", "V5", "V6")
tss <- data.frame(chr = gdf$V1, gene = gdf$V4, stringsAsFactors = FALSE)
tss$tss <-  ifelse(gdf$V5 == "+", gdf$V3, gdf$V2)
tss$start <- ifelse(tss$tss - 50000 > 0, tss$tss - 50000, 0)
tss$stop <- tss$tss + 50000
tss_idx <- makeGRangesFromDataFrame(tss, keep.extra.columns = TRUE)

adf <- peak.info
adf$site_name <- NULL 
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
 #                                                         plot.caption = element_text(vjust = 1)) +
 # labs(title = "Histogram of peaks per gene",  x = "Peaks / gene", y="Frequency") + pretty_plot()

dist <- abs(mcols(tss_idx)$tss[subjectHits(ov)] - start(atacgranges)[queryHits(ov)])
exp_dist_model <- exp(-1*dist/5000)


m <- Matrix::sparseMatrix(i = c(queryHits(ov), length(atacgranges)),
                          j = c(subjectHits(ov), length(tss_idx)),
                          x = c(exp_dist_model,0))
colnames(m) <- gdf$V4 # gene name
m <- m[,which(Matrix::colSums(m) != 0)]

counts <- data.matrix(prova_10x)
fm_genescoring <- data.matrix(t(m) %*% prova_10x)
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
cds_gs = cluster_cells(cds_gs, resolution=0.95e-3)

#class <- as.data.frame(input_cds@clusters@listData[["UMAP"]][["clusters"]])
class3 <- as.data.frame(cds_gs@clusters@listData[["UMAP"]][["clusters"]])
colnames(class) <- "CLASS"
colnames(class3) <- "CLASS"
cds_gs@colData@listData[["ATAC"]] <- class$CLASS

#plot_cells(cds_gs, color_cells_by = "ATAC", label_groups_by_cluster = FALSE)

#plot_cells(cds_gs, genes = "Rorb")

ARI(class$CLASS, class3$CLASS)
AMI(class$CLASS, class3$CLASS)




######## SIgnac #########

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)

mm10.1 <- mm10[1:21,]
#mm10.1[22,] <- mm10[64,]
mm10.1 <- Seqinfo(mm10.1$V1, seqlengths= mm10.1$V2)


mm10.1@genome[] <- "mm10"

brain_assay <- CreateChromatinAssay(
  counts = prova_10x,
  sep = c("_", "_"),
  genome = mm10.1,
  min.cells = 1
)

brain <- CreateSeuratObject(
  counts = brain_assay,
  assay = 'peaks',
  project = 'ATAC'
)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
#seqlevelsStyle(annotations) <- 'UCSC'
levels(annotations@seqnames@values) <- paste0("chr", annotations@seqnames@values)
genome(annotations) <- "mm10"
annotations@seqinfo@seqnames <- paste0("chr", annotations@seqinfo@seqnames)

Xkr4 <- annotations[annotations$gene_name == "Xkr4",]
Xkr41 <- as.data.frame(Xkr4@elementMetadata@listData)
Xkr42 <- as.data.frame(Xkr4@ranges)
Xkr4 <- cbind(Xkr41, Xkr42)
Xkr4ge <- gene_anno[gene_anno$symbol == "Xkr4",]

Annotation(brain) <- annotations

da_peaks <- FindMarkers(
  object = brain,
  ident.1 = "0",
  ident.2 = "1",
  min.pct = 0.05,
  test.use = 'LR'
)

ann <- c(annotations)




brain <- RunTFIDF(brain)
brain <- FindTopFeatures(brain, min.cutoff = 'q0')
brain <- RunSVD(object = brain)


brain <- RunUMAP(
  object = brain,
  reduction = 'lsi',
  dims = 2:30
)
brain <- FindNeighbors(
  object = brain,
  reduction = 'lsi',
  dims = 2:30
)
brain <- FindClusters(
  object = brain,
  algorithm = 3,
  resolution = 0.5,
  verbose = FALSE
)

#DimPlot(object = brain, label = TRUE)
#FeaturePlot(brain, "Adarb2")

gene.activities <- GeneActivity(brain)


brain[['RNA']] <- CreateAssayObject(counts = gene.activities)

DefaultAssay(brain) <- 'RNA'



brain <- FindVariableFeatures(brain, nfeatures = 3000)
brain <- NormalizeData(
  object = brain,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(brain$nCount_RNA)
)
brain <- ScaleData(brain)
brain <- RunPCA(brain, npcs = 30)
brain <- RunUMAP(brain, dims = 1:30, reduction.name = "umap.activity")
brain <- FindNeighbors(brain, dims = 1:30)
brain <- FindClusters(brain, resolution = 0.5, algorithm = 3)


#FeaturePlot(
  object = brain,
  features = c('Sst','Pvalb',"Gad2","Neurod6","Rorb","Syt6"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)
DefaultAssay(brain) <- 'peaks'
#DimPlot(brain, group.by = "peaks_snn_res.0.5")
#FeaturePlot(brain, features = "Apoe")


ARI(brain@meta.data[["peaks_snn_res.0.5"]], brain@meta.data[["RNA_snn_res.0.5"]])
AMI(brain@meta.data[["peaks_snn_res.0.5"]], brain@meta.data[["RNA_snn_res.0.5"]])













############################################
#fun
gene_ann_reg <- rtracklayer::readGFF("../DATA/IWBBIO_2022/Genomes/mus_musculus.GRCm38.Regulatory_Build.regulatory_features.20180516.gff.gz")
gene_ann_reg$seqid <- paste0("chr", gene_ann_reg$seqid)
gene_ann_reg <- gene_ann_reg[gene_ann_reg$seqid %in% mm10.1$V1,]

mm10.1 <- mm10[1:21,]
mm10.1 <- mm10
#mm10.1$V1[22:66] <- paste0(mm10.1$V1[22:66], ".1")
#mm10.1$V1[64] <- "chrMT"
#mm_copia <- mm10.1
#mm10.1[22,] <- mm10[64,]
mm10.1 <- Seqinfo(mm10.1$V1, seqlengths= mm10.1$V2)

mm10.1@genome[] <- "mm10"
listData <- gene_ann_reg$feature_type

peaks <- rownames(brain)
clos <- ClosestFeature(brain, peaks, annotation = ann)


closp <- clos[clos$distance > 0 & clos$distance < 40,]$type[]
closp[closp$distance > 0 & closp$distance < 40,]$type[] <- "prox"
closp[closp$distance >= 40,]$type[] <- "extra"

closp2 <- clos[clos$distance != 0,]
out_peaks <- as.character(closp2$query_region)


pie <- ggplot(clos, aes(x = factor(1), fill = factor(type))) + geom_bar(width=1)
pie + coord_polar(theta = "y") + theme(axis.title.x = element_blank(), axis.title.y = element_blank())


pie <- ggplot(closp2, aes(x = factor(1), fill = factor(type))) + geom_bar(width=1)
pie + coord_polar(theta = "y") + theme(axis.title.x = element_blank(), axis.title.y = element_blank())


labeled_peaks <- read.delim("../TMPResults/labeled_peaks/classifiedPeaks.csv")
labeled_peaks_multi <- read.csv("../TMPResults/labeled_peaks/classifiedPeaks_multiCols.csv", sep = "\t")
l_prova <- read.csv("../TMPResults/labeled_peaks/classifiedPeaks.csv", sep = c(","), header = FALSE)

########### first GAM with only promoter called only "gene" from gene annotation #######

labeled_peaks <- read.csv("../TMPResults/labeled_peaks/10x_v1_1.1.0_mousebrainP50/classifiedPeaks.csv", sep = "\t")

nmax <- max(stringr::str_count(labeled_peaks$ucsclabel, ",")) + 1

labeled_peaks <- separate(labeled_peaks, col = ucsclabel, sep = ",", into = paste0("RegFunc", seq_len(nmax)))



tss_peaks <- fData(input_cds)[6:7]
sum(!is.na(tss_peaks$gene))
sum(duplicated(tss_peaks$gene[!is.na(tss_peaks$gene)]))

genes_active_promoters <- tss_peaks[!is.na(tss_peaks$gene),]
duplicated_promoters <- genes_active_promoters[duplicated(genes_active_promoters$gene),]

tss_peaks$gene[!is.na(tss_peaks$gene)][duplicated(tss_peaks$gene[!is.na(tss_peaks$gene)])]

peaks_multi_info <- fData(input_cds)#[6:7]
peaks_multi_info$is.prom <- FALSE
peaks_multi_info <- cbind(peaks_multi_info,labeled_peaks_multi[7:61])

peaks_multi_info_tsscicero <- peaks_multi_info[!is.na(peaks_multi_info$gene),]

peaks_multi_info_tsscicero[peaks_multi_info_tsscicero$gene == "Sulf1",]

sum(peaks_multi_info_tsscicero$X0 == "")

##### check for promoters peaks######
ppp <- peaks_multi_info
ppp <- as.data.frame(peaks_multi_info)
ppp[ppp == "prom"] <- TRUE

peaks_multi_info$is.prom <- apply(ppp, 1, any)
peaks_multi_info[is.na(peaks_multi_info$is.prom),]$is.prom <- FALSE





#saveRDS(peaks_multi_info, "C:/Users/loren/IWBBIO/conns/peaks_multi_info")
#peaks_multi_info[row.names(peaks_multi_info)]$is.prom <- is.element("prom", unlist(peaks_multi_info[row.names(peaks_multi_info),]))

peaks_prom <- peaks_multi_info[peaks_multi_info$is.prom == "TRUE",]
different_prom <- peaks_prom[is.na(peaks_prom$gene),]

prom_list <- rownames(peaks_prom)
non_prom_list <- rownames(peaks_multi_info)
non_prom_list <- setdiff(non_prom_list, prom_list)

non_prom_peaks <- peaks_multi_info[non_prom_list,]


#gene_closest_to_prom <- ClosestFeature(brain, gsub("_", "-", prom_list))
#gane_closest_to_non_prom <- ClosestFeature(brain, gsub("_", "-", non_prom_list))

sub <- gene_anno[gene_anno$type == "gene",]
sub <- sub[sub$seqid != "chrM",]
gene_anno_GRanges <- makeGRangesFromDataFrame(sub, seqinfo = mm10.1, seqnames.field = "seqid", keep.extra.columns = TRUE)
Annotation(brain) <- gene_anno_GRanges

#gene_closest_to_prom[duplicated(gene_closest_to_prom$symbol),]
near_gene_closest_to_prom <- gene_closest_to_prom[gene_closest_to_prom$distance <= 1000,]
near_prom_list <- gsub("-","_",near_gene_closest_to_prom$query_region)

peaks_prom <-  peaks_prom[near_prom_list,]
peaks_prom$gene_anno <- near_gene_closest_to_prom$symbol
  

promoter_peak_table <- peaks_prom
promoter_peak_table$site_name <- rownames(promoter_peak_table)
promoter_gene_mat <-
  Matrix::sparseMatrix(j=as.numeric(factor(promoter_peak_table$site_name)),
                       i=as.numeric(factor(promoter_peak_table$gene_anno)),
                       x=1)

accessibility_mat <- exprs(input_cds)
promoter_activity_scores <- accessibility_mat[near_prom_list,, drop=FALSE]

colnames(promoter_gene_mat) = levels(factor(promoter_peak_table$site_name))
row.names(promoter_gene_mat) = levels(factor(promoter_peak_table$gene_anno))
promoter_gene_mat <- promoter_gene_mat[,row.names(promoter_activity_scores)]

#Znrf1 <- genes_active_promoters[genes_active_promoters$gene == "Znrf1",]

first_gene_matrix <- promoter_gene_mat %*% promoter_activity_scores




first_cell <- colnames(first_gene_matrix)
#cicero_cell <- read.table(cicero_cell)
lenght1 <- length(first_cell)

first_gene <- row.names(first_gene_matrix)
lenght2 <- length(first_gene)
c_c <- matrix(first_cell, nrow = lenght1, dimnames = list(first_cell,c("Cells")))
c_g<- matrix(first_gene, nrow = lenght2, dimnames = list(first_gene,c("gene_short_name")))


## processing GAM with Cicero
cds_first <-  suppressWarnings(new_cell_data_set(first_gene_matrix, cell_metadata = c_c, gene_metadata = c_g))

cds_first <- detect_genes(cds_first)
cds_first <- estimate_size_factors(cds_first)
cds_first <- preprocess_cds(cds_first, method = "PCA")

cds_first <- reduce_dimension(cds_first, reduction_method = 'UMAP', 
                               preprocess_method = "PCA")
cds_first = cluster_cells(cds_first, resolution=0.95e-3)

#plot_cells(cds_first)

cds_first@colData@listData[["ATAC"]] <- class$CLASS

#plot_cells(cds_first, color_cells_by = "ATAC", label_groups_by_cluster = FALSE)

class_first <- as.data.frame(cds_first@clusters@listData[["UMAP"]][["clusters"]])
colnames(class_first) <- "CLASS"

ARI(class$CLASS, class_first$CLASS)
AMI(class$CLASS, class_first$CLASS)



###### using refseq gene annotation ##########

refseq_gene_anno <-  rtracklayer::readGFF("../DATA/IWBBIO_2022/Genomes/refseq/GCF_000001635.26_GRCm38.p6_genomic.gtf.gz")
chr2acc <- read.csv("../DATA/IWBBIO_2022/Genomes/refseq/chr2acc.txt", sep = "\t")

refseq_peaks_prom <- peaks_multi_info[peaks_multi_info$is.prom == "TRUE",]
refseq_prom_list <- rownames(refseq_peaks_prom)

refseq_gene_anno <- refseq_gene_anno[refseq_gene_anno$seqid %in% chr2acc$Accession.version,]
refseq_gene_anno$seqid <- as.factor(as.character(refseq_gene_anno$seqid)) 
levels(refseq_gene_anno$seqid) <- chr2acc$X.Chromosome
refseq_gene_anno$seqid <- paste0("chr", refseq_gene_anno$seqid)
refseq_gene_anno <- refseq_gene_anno[refseq_gene_anno$type == "gene",]
refseq_gene_anno <- refseq_gene_anno[refseq_gene_anno$gene_biotype %in% c("protein_coding","lncRNA"),]
refseq_GRanges <- makeGRangesFromDataFrame(refseq_gene_anno, seqinfo = mm10.1, seqnames.field = "seqid", keep.extra.columns = TRUE)

refseq_gene_closest_to_prom <- ClosestFeature(brain, gsub("_", "-", refseq_prom_list), annotation = refseq_GRanges)
refseq_gene_closest_to_prom <- refseq_gene_closest_to_prom[refseq_gene_closest_to_prom$distance <= 1000,]

refseq_near_prom_list <- gsub("-","_",refseq_gene_closest_to_prom$query_region)


refseq_peaks_prom <-  refseq_peaks_prom[refseq_near_prom_list,]
refseq_peaks_prom$gene_anno <- refseq_gene_closest_to_prom$gene
refseq_peaks_prom$site_name <- rownames(refseq_peaks_prom)

refseq_promoter_peak_table <- refseq_peaks_prom[,c("gene_anno","site_name")]
#refseq_promoter_peak_table$site_name <- rownames(refseq_promoter_peak_table)
refseq_promoter_gene_mat <-
  Matrix::sparseMatrix(j=as.numeric(factor(refseq_promoter_peak_table$site_name)),
                       i=as.numeric(factor(refseq_promoter_peak_table$gene_anno)),
                       x=1)

refseq_accessibility_mat <- exprs(input_cds)
refseq_promoter_activity_scores <- refseq_accessibility_mat[refseq_near_prom_list,, drop=FALSE]

colnames(refseq_promoter_gene_mat) = levels(factor(refseq_promoter_peak_table$site_name))
row.names(refseq_promoter_gene_mat) = levels(factor(refseq_promoter_peak_table$gene_anno))
refseq_promoter_gene_mat <- refseq_promoter_gene_mat[,row.names(refseq_promoter_activity_scores)]

#Znrf1 <- genes_active_promoters[genes_active_promoters$gene == "Znrf1",]

refseq_first_gene_matrix <- refseq_promoter_gene_mat %*% refseq_promoter_activity_scores


first_cell <- colnames(refseq_first_gene_matrix)
#cicero_cell <- read.table(cicero_cell)
lenght1 <- length(first_cell)

first_gene <- row.names(refseq_first_gene_matrix)
lenght2 <- length(first_gene)
c_c <- matrix(first_cell, nrow = lenght1, dimnames = list(first_cell,c("Cells")))
c_g<- matrix(first_gene, nrow = lenght2, dimnames = list(first_gene,c("gene_short_name")))


# processing GAM with Cicero
refseq_cds_first <-  suppressWarnings(new_cell_data_set(refseq_first_gene_matrix, cell_metadata = c_c, gene_metadata = c_g))

refseq_cds_first <- detect_genes(refseq_cds_first)
refseq_cds_first <- estimate_size_factors(refseq_cds_first)
refseq_cds_first <- preprocess_cds(refseq_cds_first, method = "PCA")

refseq_cds_first <- reduce_dimension(refseq_cds_first, reduction_method = 'UMAP', 
                              preprocess_method = "PCA")
refseq_cds_first = cluster_cells(refseq_cds_first, resolution=0.95e-3)

#plot_cells(refseq_cds_first)

refseq_cds_first@colData@listData[["ATAC"]] <- class$CLASS

#plot_cells(refseq_cds_first, color_cells_by = "ATAC", label_groups_by_cluster = FALSE)

refseq_class_first <- as.data.frame(refseq_cds_first@clusters@listData[["UMAP"]][["clusters"]])
colnames(refseq_class_first) <- "CLASS"

ARI(class$CLASS, refseq_class_first$CLASS)
AMI(class$CLASS, refseq_class_first$CLASS)

confusion_sub <- choose_cells(cds_first)
confusion_sub <- detect_genes(confusion_sub)
confusion_sub <- preprocess_cds(confusion_sub, method = "PCA")
confusion_sub <- reduce_dimension(confusion_sub, reduction_method = 'UMAP', 
                                  preprocess_method = "PCA")
confusion_sub <- cluster_cells(confusion_sub, resolution=0.95e-3)
#plot_cells(confusion_sub, color_cells_by = "ATAC")

############### using only transcript  ##########

transcript_peaks_prom <- peaks_multi_info[peaks_multi_info$is.prom == "TRUE",]
transcript_prom_list <- rownames(transcript_peaks_prom)


transcript <- gene_anno[gene_anno$type == "transcript",]
transcript <- transcript[transcript$seqid != "chrM",]
transcript <- transcript[transcript$gene_type %in% c("protein_coding","lincRNA"),]
transcript <- transcript[transcript$transcript_type != "nonsense_mediated_decay",]
#refseq_gene_anno <- refseq_gene_anno[refseq_gene_anno$seqid %in% chr2acc$Accession.version,]
#refseq_gene_anno$seqid <- as.factor(as.character(refseq_gene_anno$seqid)) 
#levels(refseq_gene_anno$seqid) <- chr2acc$X.Chromosome
#refseq_gene_anno$seqid <- paste0("chr", refseq_gene_anno$seqid)
#refseq_gene_anno <- refseq_gene_anno[refseq_gene_anno$type == "gene",]

transcript_GRanges <- makeGRangesFromDataFrame(transcript, seqinfo = mm10.1, seqnames.field = "seqid", keep.extra.columns = TRUE)

transcript_gene_closest_to_prom <- ClosestFeature(brain, gsub("_", "-", refseq_prom_list), annotation = transcript_GRanges)
transcript_gene_closest_to_prom <- transcript_gene_closest_to_prom[transcript_gene_closest_to_prom$distance <= 1000,]

transcript_near_prom_list <- gsub("-","_",transcript_gene_closest_to_prom$query_region)


transcript_peaks_prom <-  transcript_peaks_prom[transcript_near_prom_list,]
transcript_peaks_prom$gene_anno <- transcript_gene_closest_to_prom$symbol


transcript_promoter_peak_table <- transcript_peaks_prom
transcript_promoter_peak_table$site_name <- rownames(transcript_promoter_peak_table)
transcript_promoter_gene_mat <-
  Matrix::sparseMatrix(j=as.numeric(factor(transcript_promoter_peak_table$site_name)),
                       i=as.numeric(factor(transcript_promoter_peak_table$gene_anno)),
                       x=1)

transcript_accessibility_mat <- exprs(input_cds)
transcript_promoter_activity_scores <- transcript_accessibility_mat[transcript_near_prom_list,, drop=FALSE]

colnames(transcript_promoter_gene_mat) = levels(factor(transcript_promoter_peak_table$site_name))
row.names(transcript_promoter_gene_mat) = levels(factor(transcript_promoter_peak_table$gene_anno))
transcript_promoter_gene_mat <- transcript_promoter_gene_mat[,row.names(transcript_promoter_activity_scores)]

#Znrf1 <- genes_active_promoters[genes_active_promoters$gene == "Znrf1",]

transcript_first_gene_matrix <- transcript_promoter_gene_mat %*% transcript_promoter_activity_scores


first_cell <- colnames(transcript_first_gene_matrix)
#cicero_cell <- read.table(cicero_cell)
lenght1 <- length(first_cell)

first_gene <- row.names(transcript_first_gene_matrix)
lenght2 <- length(first_gene)
c_c <- matrix(first_cell, nrow = lenght1, dimnames = list(first_cell,c("Cells")))
c_g<- matrix(first_gene, nrow = lenght2, dimnames = list(first_gene,c("gene_short_name")))


## processing GAM with Cicero
transcript_cds_first <-  suppressWarnings(new_cell_data_set(transcript_first_gene_matrix, cell_metadata = c_c, gene_metadata = c_g))

transcript_cds_first <- detect_genes(transcript_cds_first)
transcript_cds_first <- estimate_size_factors(transcript_cds_first)
transcript_cds_first <- preprocess_cds(transcript_cds_first, method = "PCA")

transcript_cds_first <- reduce_dimension(transcript_cds_first, reduction_method = 'UMAP', 
                                     preprocess_method = "PCA")
transcript_cds_first = cluster_cells(transcript_cds_first, resolution=0.95e-3)

#plot_cells(transcript_cds_first)

transcript_cds_first@colData@listData[["ATAC"]] <- class$CLASS

#plot_cells(transcript_cds_first, color_cells_by = "ATAC", label_groups_by_cluster = FALSE)

transcript_class_first <- as.data.frame(transcript_cds_first@clusters@listData[["UMAP"]][["clusters"]])
colnames(transcript_class_first) <- "CLASS"

ARI(class$CLASS, transcript_class_first$CLASS)
AMI(class$CLASS, transcript_class_first$CLASS)


########### adding extra intergenetic peaks to GAM ############

peaks_multi_info
peaks_prom <- peaks_multi_info[peaks_multi_info$is.prom == "TRUE",]

#different_prom <- peaks_prom[is.na(peaks_prom$gene),]

prom_list <- rownames(peaks_prom)
#non_prom_list <- rownames(peaks_multi_info)
#non_prom_list <- setdiff(non_prom_list, prom_list)
non_prom_peaks <- peaks_multi_info[peaks_multi_info$is.prom != "TRUE",]
non_prom_list <- rownames(non_prom_peaks)

intragenetic_non_prom_peaks <- ClosestFeature(brain, gsub("_", "-", non_prom_list), refseq_GRanges)
intragenetic_non_prom_peaks <- intragenetic_non_prom_peaks[intragenetic_non_prom_peaks$distance == 0,]
intragenetic_non_prom_peaks_list <- gsub("-","_",intragenetic_non_prom_peaks$query_region)


intragenetic_peaks <-  non_prom_peaks[intragenetic_non_prom_peaks_list,]
intragenetic_peaks$gene_anno <- intragenetic_non_prom_peaks$gene
intragenetic_peaks$site_name <- rownames(intragenetic_peaks)

intragenetic_peak_table <- intragenetic_peaks[,c("gene_anno","site_name")]
refseq_promoter_peak_table

inta_gene_name <- levels(factor(intragenetic_peak_table$gene_anno))
prom_gene_name <- levels(factor(refseq_promoter_peak_table$gene_anno))

diff_gene_name <- setdiff(inta_gene_name,prom_gene_name)

peak_table <- rbind(refseq_promoter_peak_table, intragenetic_peak_table)
allgene_peaks_list <- rownames(peak_table)
#intragenetic_peak_table$site_name <- rownames(intragenetic_peaks)
refseq_allgene_gene_mat <-
  Matrix::sparseMatrix(j=as.numeric(factor(peak_table$site_name)),
                       i=as.numeric(factor(peak_table$gene_anno)),
                       x=1)

refseq_allgene_accessibility_mat <- exprs(input_cds)
refseq_allgene_activity_scores <- refseq_allgene_accessibility_mat[allgene_peaks_list,, drop=FALSE]

colnames(refseq_allgene_gene_mat) = levels(factor(peak_table$site_name))
row.names(refseq_allgene_gene_mat) = levels(factor(peak_table$gene_anno))
refseq_allgene_gene_mat <- refseq_allgene_gene_mat[,row.names(refseq_allgene_activity_scores)]

refseq_allgene_gene_matrix <- refseq_allgene_gene_mat %*% refseq_allgene_activity_scores


first_cell <- colnames(refseq_allgene_gene_matrix)
lenght1 <- length(first_cell)
first_gene <- row.names(refseq_allgene_gene_matrix)
lenght2 <- length(first_gene)
c_c <- matrix(first_cell, nrow = lenght1, dimnames = list(first_cell,c("Cells")))
c_g<- matrix(first_gene, nrow = lenght2, dimnames = list(first_gene,c("gene_short_name")))


# processing GAM with Cicero
refseq_allgene_cds_first <-  suppressWarnings(new_cell_data_set(refseq_allgene_gene_matrix, cell_metadata = c_c, gene_metadata = c_g))

refseq_allgene_cds_first <- detect_genes(refseq_allgene_cds_first)
refseq_allgene_cds_first <- estimate_size_factors(refseq_allgene_cds_first)
refseq_allgene_cds_first <- preprocess_cds(refseq_allgene_cds_first, method = "PCA")

refseq_allgene_cds_first <- reduce_dimension(refseq_allgene_cds_first, reduction_method = 'UMAP', 
                                     preprocess_method = "PCA")
refseq_allgene_cds_first = cluster_cells(refseq_allgene_cds_first, resolution=0.95e-3)

#plot_cells(refseq_allgene_cds_first)

refseq_allgene_cds_first@colData@listData[["ATAC"]] <- class$CLASS

#plot_cells(refseq_allgene_cds_first, color_cells_by = "ATAC", label_groups_by_cluster = FALSE)

refseq_allgene_class_first <- as.data.frame(refseq_allgene_cds_first@clusters@listData[["UMAP"]][["clusters"]])
colnames(refseq_allgene_class_first) <- "CLASS"

ARI(class$CLASS, refseq_allgene_class_first$CLASS)
AMI(class$CLASS, refseq_allgene_class_first$CLASS)


######### consider only genes with active promoters #####

inta_gene_name <- levels(factor(intragenetic_peak_table$gene_anno))
prom_gene_name <- levels(factor(refseq_promoter_peak_table$gene_anno))

diff_gene_name <- setdiff(inta_gene_name,prom_gene_name)
#no active promoter genes identified
length(diff_gene_name)

#start with peak_table

refseq_allgene_gene_mat <-
  Matrix::sparseMatrix(j=as.numeric(factor(peak_table$site_name)),
                       i=as.numeric(factor(peak_table$gene_anno)),
                       x=1)
only_prom <- refseq_allgene_gene_mat[prom_gene_name,]


only_prom_accessibility_mat <- exprs(input_cds)
only_prom_activity_scores <- only_prom_accessibility_mat[allgene_peaks_list,, drop=FALSE]

colnames(only_prom) = levels(factor(peak_table$site_name))
row.names(only_prom) = levels(factor(refseq_promoter_peak_table$gene_anno))
only_prom <- only_prom[,row.names(only_prom_activity_scores)]

only_prom_gene_matrix <- only_prom %*% only_prom_activity_scores



first_cell <- colnames(only_prom_gene_matrix)
lenght1 <- length(first_cell)
first_gene <- row.names(only_prom_gene_matrix)
lenght2 <- length(first_gene)
c_c <- matrix(first_cell, nrow = lenght1, dimnames = list(first_cell,c("Cells")))
c_g<- matrix(first_gene, nrow = lenght2, dimnames = list(first_gene,c("gene_short_name")))


# processing GAM with Cicero
only_prom_cds_first <-  suppressWarnings(new_cell_data_set(only_prom_gene_matrix, cell_metadata = c_c, gene_metadata = c_g))

only_prom_cds_first <- detect_genes(only_prom_cds_first)
only_prom_cds_first <- estimate_size_factors(only_prom_cds_first)
only_prom_cds_first <- preprocess_cds(only_prom_cds_first, method = "PCA")

only_prom_cds_first <- reduce_dimension(only_prom_cds_first, reduction_method = 'UMAP', 
                                             preprocess_method = "PCA")
only_prom_cds_first = cluster_cells(only_prom_cds_first, resolution=0.95e-3)

#plot_cells(only_prom_cds_first)

only_prom_cds_first@colData@listData[["ATAC"]] <- class$CLASS

#plot_cells(only_prom_cds_first, color_cells_by = "ATAC", label_groups_by_cluster = FALSE)

only_prom_class_first <- as.data.frame(only_prom_cds_first@clusters@listData[["UMAP"]][["clusters"]])
colnames(only_prom_class_first) <- "CLASS"

ARI(class$CLASS, only_prom_class_first$CLASS)
AMI(class$CLASS, only_prom_class_first$CLASS)

GA_only_prom <- top_markers(only_prom_cds_first)



################# normalization attempt ############
#cicero normalization 
unnorm_first <- first_gene_matrix[!Matrix::rowSums(first_gene_matrix) == 0, 
                       !Matrix::colSums(first_gene_matrix) == 0]
num_genes <- pData(input_cds)$num_genes_expressed
names(num_genes) <- row.names(pData(input_cds))

norm_first_gene_matrix <- normalize_gene_activities(unnorm_first, num_genes)

norm_first_cell <- colnames(norm_first_gene_matrix)
  #cicero_cell <- read.table(cicero_cell)
lenght1 <- length(norm_first_cell)

norm_first_gene <- row.names(norm_first_gene_matrix)
lenght2 <- length(norm_first_gene)
c_c <- matrix(norm_first_cell, nrow = lenght1, dimnames = list(norm_first_cell,c("Cells")))
c_g<- matrix(norm_first_gene, nrow = lenght2, dimnames = list(norm_first_gene,c("gene_short_name")))


## processing GAM with Cicero
norm_cds_first <-  suppressWarnings(new_cell_data_set(norm_first_gene_matrix, cell_metadata = c_c, gene_metadata = c_g))

norm_cds_first <- detect_genes(norm_cds_first)
norm_cds_first <- estimate_size_factors(norm_cds_first)
norm_cds_first <- preprocess_cds(norm_cds_first, method = "PCA")

norm_cds_first <- reduce_dimension(norm_cds_first, reduction_method = 'UMAP', 
                              preprocess_method = "PCA")
norm_cds_first = cluster_cells(norm_cds_first, resolution=0.95e-3)

#plot_cells(norm_cds_first)

norm_cds_first@colData@listData[["ATAC"]] <- class$CLASS

#plot_cells(norm_cds_first, color_cells_by = "ATAC", label_groups_by_cluster = FALSE)

class_first_norm <- as.data.frame(norm_cds_first@clusters@listData[["UMAP"]][["clusters"]])
colnames(class_first_norm) <- "CLASS"

ARI(class$CLASS, class_first_norm$CLASS)
AMI(class$CLASS, class_first_norm$CLASS)

################# refseqFuncElements ############

library(dplyr)
library(tidyr)
library(stringr)

refseq_funcelement <- read.csv("../DATA/IWBBIO_2022/Genomes/refseq/refseq FuncElement/refSeqFuncElems_soTerm_classifiedPeaks.csv", sep = "\t")
refseq_funcelement_names <- read.csv("../DATA/IWBBIO_2022/Genomes/refseq/refseq FuncElement/refSeqFuncElems_name_classifiedPeaks.csv", sep = "\t")

nmax <- max(stringr::str_count(refseq_funcelement$refSeqFuncElems_soTerm, "\t")) + 1

refseq_funcelement <- separate(refseq_funcelement, col = refSeqFuncElems_soTerm, sep = "\t", into = paste0("FuncElem", seq_len(nmax)))


nmax <- max(stringr::str_count(refseq_funcelement_names$refSeqFuncElems_soTerm, "\t")) + 1

refseq_funcelement_names <- separate(refseq_funcelement_names, col = refSeqFuncElems_soTerm, sep = "\t", into = paste0("FuncElem", seq_len(nmax)))

#refseq_funcelement[refseq_funcelement$FuncElem1 != "",]

ppp <- as.data.frame(refseq_funcelement)
ppp[ppp == "enhancer"] <- TRUE

peaks_multi_info$is.enh <- FALSE
peaks_multi_info$is.enh <- apply(ppp, 1, any)

peaks_multi_info[is.na(peaks_multi_info$is.enh),]$is.enh <- FALSE
enh_peaks_list <- rownames(peaks_multi_info[peaks_multi_info$is.enh == TRUE,])

length(setdiff(enh_peaks_list, prom_list))

ppp <- as.data.frame(peaks_multi_info)

############## adding connection ########### 

dist_thresh=300000

accessibility_mat <- exprs(input_cds)
promoter_peak_table <- fData(input_cds)
promoter_peak_table <- refseq_peaks_prom[, c("site_name", "gene_anno")]


promoter_peak_table$site_name <- as.character(row.names(promoter_peak_table))
promoter_peak_table <- promoter_peak_table[!is.na(promoter_peak_table$gene_anno),]
promoter_peak_table <- promoter_peak_table[,c("site_name", "gene_anno")]
promoter_peak_table$gene_anno <- as.character(promoter_peak_table$gene_anno)

colnames(promoter_peak_table) <- c("peak", "gene")
  
split_peak_names <- function(inp) {
  out <- stringr::str_split_fixed(stringi::stri_reverse(inp), 
                                  ":|-|_", 3)
  out[,1] <- stringi::stri_reverse(out[,1])
  out[,2] <- stringi::stri_reverse(out[,2])
  out[,3] <- stringi::stri_reverse(out[,3])
  out[,c(3,2,1), drop=FALSE]
}

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


{
nonneg_cons <-
  conns[(conns$Peak1 %in%
                      promoter_peak_table$peak |
           conns$Peak2 %in%
                      promoter_peak_table$peak) &
          conns$coaccess >= coaccess &
          conns$dist < dist_thresh,]
nonneg_cons <- nonneg_cons[,c("Peak1", "Peak2", "coaccess")]
nonneg_cons <- nonneg_cons[!duplicated(nonneg_cons),]

nonneg_cons$Peak1 <- as.character(nonneg_cons$Peak1)
nonneg_cons$Peak2 <- as.character(nonneg_cons$Peak2)

nonneg_cons <- rbind(nonneg_cons,
                     data.frame(Peak1=unique(promoter_peak_table$peak),
                                Peak2=unique(promoter_peak_table$peak),
                                coaccess=0))
}

prom_enh <- conns[(conns$Peak1 %in%
         refseq_near_prom_list &
         conns$Peak2 %in%
         enh_peaks_list),]

prom_enh <- conns[(conns$Peak1 %in%
                     promoter_peak_table$peak &
                     conns$Peak2 %in%
                     enh_peaks_list) | (conns$Peak2 %in%
                                          promoter_peak_table$peak &
                                          conns$Peak1 %in%
                                          enh_peaks_list),]

prom_enh <- prom_enh[!duplicated(prom_enh),]
prom_enh <- prom_enh[prom_enh$coaccess >= 0.12 & prom_enh$dist <= dist_thresh,]

prom_enh <- prom_enh[,c("Peak1", "Peak2", "coaccess")]
prom_enh <- prom_enh[!duplicated(prom_enh),]

prom_enh$Peak1 <- as.character(prom_enh$Peak1)
prom_enh$Peak2 <- as.character(prom_enh$Peak2)


prom_enh <- rbind(prom_enh,
                     data.frame(Peak1=unique(promoter_peak_table$peak),
                                Peak2=unique(promoter_peak_table$peak),
                                coaccess=0))


prom_enh_connectivity <- make_sparse_matrix(prom_enh, x.name = "coaccess")
#distal_connectivity_matrix <- make_sparse_matrix(nonneg_cons, x.name="coaccess")




promoter_conn_matrix <-
  prom_enh_connectivity[unique(promoter_peak_table$peak),]

# Get list of promoter and distal sites in accessibility mat
promoter_safe_sites <- intersect(rownames(promoter_conn_matrix),
                                 row.names(accessibility_mat))
distal_safe_sites <- intersect(colnames(promoter_conn_matrix),
                               row.names(accessibility_mat))
distal_safe_sites <- setdiff(distal_safe_sites, promoter_safe_sites)

# Get accessibility info for promoters
promoter_access_mat_in_cicero_map <- accessibility_mat[promoter_safe_sites,, drop=FALSE]

# Get accessibility for distal sites
distal_activity_scores <- accessibility_mat[distal_safe_sites,, drop=FALSE]

# Scale connectivity matrix by site_weights
scaled_site_weights <- site_weights[distal_safe_sites,distal_safe_sites, drop=FALSE]
total_linked_site_weights <- promoter_conn_matrix[,distal_safe_sites, drop=FALSE] %*%
  scaled_site_weights
total_linked_site_weights <- 1/Matrix::rowSums(total_linked_site_weights,
                                               na.rm=TRUE)
total_linked_site_weights[is.finite(total_linked_site_weights) == FALSE] <- 0
total_linked_site_weights[is.na(total_linked_site_weights)] <- 0
total_linked_site_weights[is.nan(total_linked_site_weights)] <- 0
total_linked_site_weights <- Matrix::Diagonal(x=total_linked_site_weights)
scaled_site_weights <- total_linked_site_weights %*%
  promoter_conn_matrix[,distal_safe_sites, drop=FALSE] %*%promoter_activity_scores
  scaled_site_weights
scaled_site_weights@x[scaled_site_weights@x > 1] <- 1

# Multiply distal accessibility by site weights
distal_activity_scores <- scaled_site_weights %*% distal_activity_scores

distal_activity_scores <-
  distal_activity_scores[row.names(promoter_access_mat_in_cicero_map),, drop=FALSE]

# Sum distal and promoter scores
promoter_activity_scores <- distal_activity_scores +  
  promoter_access_mat_in_cicero_map



# Make and populate final matrix
promoter_gene_mat <-
  Matrix::sparseMatrix(j=as.numeric(factor(promoter_peak_table$peak)),
                       i=as.numeric(factor(promoter_peak_table$gene)),
                       x=1)
colnames(promoter_gene_mat) = levels(factor(promoter_peak_table$peak))
row.names(promoter_gene_mat) = levels(factor(promoter_peak_table$gene))
promoter_gene_mat <- promoter_gene_mat[,row.names(promoter_activity_scores)]
gene_activity_scores <- promoter_gene_mat %*% promoter_activity_scores


first_cell <- colnames(gene_activity_scores)
#cicero_cell <- read.table(cicero_cell)
lenght1 <- length(first_cell)

first_gene <- row.names(gene_activity_scores)
lenght2 <- length(first_gene)
c_c <- matrix(first_cell, nrow = lenght1, dimnames = list(first_cell,c("Cells")))
c_g<- matrix(first_gene, nrow = lenght2, dimnames = list(first_gene,c("gene_short_name")))


# processing GAM with Cicero
refseq_conn_cds <-  suppressWarnings(new_cell_data_set(gene_activity_scores, cell_metadata = c_c, gene_metadata = c_g))

refseq_conn_cds <- detect_genes(refseq_conn_cds)
refseq_conn_cds <- estimate_size_factors(refseq_conn_cds)
refseq_conn_cds <- preprocess_cds(refseq_conn_cds, method = "PCA")

refseq_conn_cds <- reduce_dimension(refseq_conn_cds, reduction_method = 'UMAP', 
                                     preprocess_method = "PCA")
refseq_conn_cds = cluster_cells(refseq_conn_cds, resolution=0.95e-3)

#plot_cells(refseq_conn_cds)

refseq_conn_cds@colData@listData[["ATAC"]] <- class$CLASS

#plot_cells(refseq_conn_cds, color_cells_by = "ATAC", label_groups_by_cluster = FALSE)

refseq_conn_class_first <- as.data.frame(refseq_conn_cds@clusters@listData[["UMAP"]][["clusters"]])
colnames(refseq_conn_class_first) <- "CLASS"

ARI(class$CLASS, refseq_conn_class_first$CLASS)
AMI(class$CLASS, refseq_conn_class_first$CLASS)

############# CCRE enhd conns #######

ppp <- peaks_multi_info
ppp <- as.data.frame(peaks_multi_info)
ppp[ppp == "enhD"] <- TRUE

peaks_multi_info$is.enhD <- FALSE
peaks_multi_info$is.enhD <- apply(ppp[4:58],1, any)
peaks_multi_info[is.na(peaks_multi_info$is.enhD),]$is.enhD <- FALSE

ppp <- peaks_multi_info[peaks_multi_info$is.prom == "TRUE",]

enhD_peaks_list <- rownames(peaks_multi_info[peaks_multi_info$is.enhD == TRUE,])
enhD_peaks_list <- setdiff(enhD_peaks_list, prom_list)


promoter_peak_table <- refseq_peaks_prom[, c("site_name", "gene_anno")]


promoter_peak_table$site_name <- as.character(row.names(promoter_peak_table))
promoter_peak_table <- promoter_peak_table[!is.na(promoter_peak_table$gene_anno),]
promoter_peak_table <- promoter_peak_table[,c("site_name", "gene_anno")]
promoter_peak_table$gene_anno <- as.character(promoter_peak_table$gene_anno)

colnames(promoter_peak_table) <- c("peak", "gene")


prom_enhD <- conns[(conns$Peak1 %in%
                     promoter_peak_table$peak &
                     conns$Peak2 %in%
                     enhD_peaks_list) | (conns$Peak2 %in%
                                          promoter_peak_table$peak &
                                          conns$Peak1 %in%
                                          enhD_peaks_list),]
prom_enhD <- prom_enhD[prom_enhD$coaccess >= 0.12 & prom_enhD$dist <= dist_thresh,]
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






####### adding enhD to refseq prom score

GAM.1 <- promoter_gene_mat2 %*% distal_activity_scores2
GAM.2 <- only_prom_gene_matrix

GAM_final <- GAM.1 +GAM.2



first_cell <- colnames(GAM_final)
#cicero_cell <- read.table(cicero_cell)
lenght1 <- length(first_cell)

first_gene <- row.names(GAM_final)
lenght2 <- length(first_gene)
c_c <- matrix(first_cell, nrow = lenght1, dimnames = list(first_cell,c("Cells")))
c_g<- matrix(first_gene, nrow = lenght2, dimnames = list(first_gene,c("gene_short_name")))


# processing GAM with Cicero
GAM_final_cds <-  suppressWarnings(new_cell_data_set(GAM_final, cell_metadata = c_c, gene_metadata = c_g))

GAM_final_cds <- detect_genes(GAM_final_cds)
GAM_final_cds <- estimate_size_factors(GAM_final_cds)
GAM_final_cds <- preprocess_cds(GAM_final_cds, method = "PCA")

GAM_final_cds <- reduce_dimension(GAM_final_cds, reduction_method = 'UMAP', 
                                    preprocess_method = "PCA")
GAM_final_cds = cluster_cells(GAM_final_cds, resolution=0.95e-3)

#plot_cells(GAM_final_cds)

final_GAM_cds@colData@listData[["ATAC"]] <- class$CLASS

#plot_cells(GAM_final_cds, color_cells_by = "ATAC", label_groups_by_cluster = FALSE)

GAM_final_class <- as.data.frame(GAM_final_cds@clusters@listData[["UMAP"]][["clusters"]])
colnames(GAM_final_class) <- "CLASS"

ARI(class$CLASS, GAM_final_class$CLASS)
AMI(class$CLASS, GAM_final_class$CLASS)



g8 <- choose_cells(GAM_final_cds)

markers <- top_markers(g8, group_cells_by = "ATAC")


#plot_cells(g8,  color_cells_by = "ATAC", label_groups_by_cluster = FALSE)
#plot_cells(g8, genes = "Rorb")

###### norm #########
unnorm_final <- GAM_final[!Matrix::rowSums(GAM_final) == 0, 
                                  !Matrix::colSums(GAM_final) == 0]
num_genes <- pData(input_cds)$num_genes_expressed
names(num_genes) <- row.names(pData(input_cds))

norm_final_gene_matrix <- normalize_gene_activities(unnorm_final, num_genes)

norm_first_cell <- colnames(norm_final_gene_matrix)
#cicero_cell <- read.table(cicero_cell)
lenght1 <- length(norm_first_cell)

norm_first_gene <- row.names(norm_final_gene_matrix)
lenght2 <- length(norm_first_gene)
c_c <- matrix(norm_first_cell, nrow = lenght1, dimnames = list(norm_first_cell,c("Cells")))
c_g<- matrix(norm_first_gene, nrow = lenght2, dimnames = list(norm_first_gene,c("gene_short_name")))


## processing GAM with Cicero
norm_final_cds_first <-  suppressWarnings(new_cell_data_set(norm_final_gene_matrix, cell_metadata = c_c, gene_metadata = c_g))

norm_final_cds_first <- detect_genes(norm_final_cds_first)
norm_final_cds_first <- estimate_size_factors(norm_final_cds_first)
norm_final_cds_first <- preprocess_cds(norm_final_cds_first, method = "PCA")

norm_final_cds_first <- reduce_dimension(norm_final_cds_first, reduction_method = 'UMAP', 
                                   preprocess_method = "PCA")
norm_final_cds_first = cluster_cells(norm_final_cds_first, resolution=0.95e-3)

#plot_cells(norm_final_cds_first)

norm_final_cds_first@colData@listData[["ATAC"]] <- class$CLASS

#plot_cells(norm_final_cds_first, color_cells_by = "ATAC", label_groups_by_cluster = FALSE)








































































































class_final_norm <- as.data.frame(norm_final_cds_first@clusters@listData[["UMAP"]][["clusters"]])
colnames(class_final_norm) <- "CLASS"

ARI(class$CLASS, class_final_norm$CLASS)
AMI(class$CLASS, class_final_norm$CLASS)


#########
labeled_peaks_multi_sub <- labeled_peaks_multi[7:61]
cbind(labeled_peaks_multi_sub[1:55])
labeled_peaks_multi_sub <- as.vector(labeled_peaks_multi_sub)
lp <- labeled_peaks_multi[8:61]

#lp <- c(sum(lp == "enhD"), sum(lp == "enhP"), sum(lp == "CTCF"),  sum(lp == "prom"), sum(lp == "K4m3"))
#lp <- as.data.frame(lp)

concat <- labeled_peaks_multi[7]
colnames(concat) <- "V1"
lp <- lp[lp!= ""]
concat <- bind_cols(concat, lp)
for (i in colnames(lp)){
  concat <- rbind(concat, lp[lp[i]!= ""][i])
  
  
}






pie <- ggplot(labeled_peaks_multi, aes(x = factor(1), fill = factor(X0))) + geom_bar(width=1)
pie + coord_polar(theta = "y") + theme(axis.title.x = element_blank(), axis.title.y = element_blank())

pie <- ggplot(lp, aes(x = factor(1), fill = factor(X0))) + geom_bar(width=1)
pie + coord_polar(theta = "y") + theme(axis.title.x = element_blank(), axis.title.y = element_blank())



labeled_peaks[,7] <- strsplit(x = labeled_peaks[,7], split = ",")
lp <- as.vector(strsplit(x = labeled_peaks$ucsclabel, split = ","))
l <-  as.data.frame(matrix(unlist(lp), nrow=length(unlist(lp[1]))))
as.data.frame.list(lp)
########################

prova <- makeGRangesFromDataFrame(gene_ann_reg, seqinfo = mm10.1, seqnames.field = "seqid", keep.extra.columns = TRUE)


clos23 <- ClosestFeature(brain, peaks, annotation = prova2)


closp23 <- clos23[clos23$distance <= 50,]


pie <- ggplot(clos, aes(x = factor(1), fill = factor(type))) + geom_bar(width=1)
pie + coord_polar(theta = "y") + theme(axis.title.x = element_blank(), axis.title.y = element_blank())


pie <- ggplot(closp23, aes(x = factor(1), fill = factor(feature_type))) + geom_bar(width=1)
pie + coord_polar(theta = "y") + theme(axis.title.x = element_blank(), axis.title.y = element_blank())

colnames(gene_ann_reg)[3] <- "what"
colnames(gene_ann_reg)[13] <- "type"

prova2 <- bind_rows(gene_anno, gene_ann_reg, )
prova2 <- prova2[prova2$seqid != "chrM",]
prova2 <- makeGRangesFromDataFrame(prova2, seqinfo = mm10.1, seqnames.field = "seqid", keep.extra.columns = TRUE)


clos234 <- ClosestFeature(brain, peaks, annotation = prova2)


closp234 <- clos234[clos234$distance <= 50,]


pie <- ggplot(clos234, aes(x = factor(1), fill = factor(type))) + geom_bar(width=1)
pie + coord_polar(theta = "y") + theme(axis.title.x = element_blank(), axis.title.y = element_blank())


pie <- ggplot(closp234, aes(x = factor(1), fill = factor(type))) + geom_bar(width=1)
pie + coord_polar(theta = "y") + theme(axis.title.x = element_blank(), axis.title.y = element_blank())


clos_anno <- c(annotations, prova)

clos_ <- ClosestFeature(brain, peaks, annotation = clos_anno)
closp_ <- clos_[clos_$distance <= 50,]
pie <- ggplot(clos_, aes(x = factor(1), fill = factor(type))) + geom_bar(width=1)
pie + coord_polar(theta = "y") + theme(axis.title.x = element_blank(), axis.title.y = element_blank())

pie <- ggplot(closp_, aes(x = factor(1), fill = factor(type))) + geom_bar(width=1)
pie + coord_polar(theta = "y") + theme(axis.title.x = element_blank(), axis.title.y = element_blank())


extragenetic_clos_regulation <- ClosestFeature(brain, out_peaks, annotation = prova)
extragenetic_clos_regulation

pie <- ggplot(extragenetic_clos_regulation, aes(x = factor(1), fill = factor(type))) + geom_bar(width=1)
pie + coord_polar(theta = "y") + theme(axis.title.x = element_blank(), axis.title.y = element_blank())

extragenetic_clos_regulation_zd <- extragenetic_clos_regulation[extragenetic_clos_regulation$distance == 0,] 

pie <- ggplot(extragenetic_clos_regulation_zd, aes(x = factor(1), fill = factor(type))) + geom_bar(width=1)
pie + coord_polar(theta = "y") + theme(axis.title.x = element_blank(), axis.title.y = element_blank())

###########
# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

# add the gene information to the object
Annotation(brain) <- annotations



transcripts <- CollapseToLongestTranscript(ranges = annotations)
if (!is.null(x = biotypes)) {
  transcripts <- transcripts[transcripts$gene_biotype %in% biotypes]
  if (length(x = transcripts) == 0) {
    stop("No genes remaining after filtering for requested biotypes")
  }
}

biotypes = "protein_coding"
features = NULL
# filter genes if provided
if (!is.null(x = features)) {
  transcripts <- transcripts[transcripts$gene_name %in% features]
  if (length(x = transcripts) == 0) {
    stop("None of the requested genes were found in the gene annotation")
  }
}
max.width = 500000
if (!is.null(x = max.width)) {
  transcript.keep <- which(x = width(x = transcripts) < max.width)
  transcripts <- transcripts[transcript.keep]
  if (length(x = transcripts) == 0) {
    stop("No genes remaining after filtering for max.width")
  }
}

extend.upstream = 2000
extend.downstream = 0

transcripts <- Extend(
  x = transcripts,
  upstream = extend.upstream,
  downstream = extend.downstream
)

frags <- Fragments(object = brain[["peaks"]])
if (length(x = frags) == 0) {
  stop("No fragment information found for requested assay")
}

cells <- colnames(x = brain[["peaks"]])
counts <- FeatureMatrix(
  fragments = frags,
  features = transcripts,
  cells = cells,
  verbose = TRUE
)


if (!is.null(x = cells)) {
  obj.use <- c()
  for (i in seq_along(along.with = frags)) {
    if (any(cells %in% Cells(x = frags[[i]]))) {
      obj.use <- c(obj.use, i)
    }
  }
} else {
  obj.use <- seq_along(along.with = frags)
}
features = transcripts
mat.list <- sapply(
  X = obj.use,
  FUN = function(x) {
    SingleFeatureMatrix(
      fragment = frags[[x]],
      features = features,
      cells = cells,
      sep = sep,
      verbose = verbose,
      process_n = process_n
    )
  })

fragment.path <- frags[[1]]@path
frag.cells <- frags[[1]]@cells
############
fragment <- frags
SingleFeatureMatrix <- function(
  fragment,
  features,
  cells = NULL,
  process_n = 2000,
  sep = c("-", "-"),
  verbose = TRUE
) {
  fragment.path <- GetFragmentData(object = fragment, slot = "path")
  if (!is.null(cells)) {
    # only look for cells that are in the fragment file
    frag.cells <- GetFragmentData(object = fragment, slot = "cells")
    # first subset frag.cells
    cell.idx <- fmatch(
      x = names(x = frag.cells),
      table = cells,
      nomatch = 0L
    ) > 0
    cells <- frag.cells[cell.idx]
  }
  tbx <- TabixFile(file = fragment.path)
  features <- keepSeqlevels(
    x = features,
    value = intersect(
      x = seqnames(x = features),
      y = seqnamesTabix(file = tbx)
    ),
    pruning.mode = "coarse"
  )
  if (length(x = features) == 0) {
    stop("No matching chromosomes found in fragment file.")
  }
  
  feature.list <- ChunkGRanges(
    granges = features,
    nchunk = ceiling(x = length(x = features) / process_n)
  )
  if (verbose) {
    message("Extracting reads overlapping genomic regions")
  }
  if (nbrOfWorkers() > 1) {
    matrix.parts <- future_lapply(
      X = feature.list,
      FUN = PartialMatrix,
      tabix = tbx,
      cells = cells,
      sep = sep,
      future.globals = list(),
      future.scheduling = FALSE
    )
  } else {
    mylapply <- ifelse(test = verbose, yes = pblapply, no = lapply)
    matrix.parts <- mylapply(
      X = feature.list,
      FUN = PartialMatrix,
      tabix = tbx,
      cells = cells,
      sep = sep
    )
  }
  # remove any that are NULL (no fragments for any cells in the region)
  null.parts <- sapply(X = matrix.parts, FUN = is.null)
  matrix.parts <- matrix.parts[!null.parts]
  if (is.null(x = cells)) {
    all.cells <- unique(
      x = unlist(x = lapply(X = matrix.parts, FUN = colnames))
    )
    matrix.parts <- lapply(
      X = matrix.parts,
      FUN = AddMissingCells,
      cells = all.cells
    )
  }
  featmat <- do.call(what = rbind, args = matrix.parts)
  if (!is.null(x = cells)) {
    # cells supplied, rename with cell name from object rather than file
    cell.convert <- names(x = cells)
    names(x = cell.convert) <- cells
    colnames(x = featmat) <- unname(obj = cell.convert[colnames(x = featmat)])
  }
  # reorder features
  feat.str <- GRangesToString(grange = features, sep = sep)
  featmat <- featmat[feat.str, ]
  return(featmat)
}

##################
CollapseToLongestTranscript <- function(ranges) {
  range.df <- as.data.table(x = ranges)
  range.df$strand <- as.character(x = range.df$strand)
  range.df$strand <- ifelse(
    test = range.df$strand == "*",
    yes = "+",
    no = range.df$strand
  )
  collapsed <- range.df[
    , .(unique(seqnames),
        min(start),
        max(end),
        strand[[1]],
        gene_biotype[[1]],
        gene_name[[1]]),
    "gene_id"
  ]
  colnames(x = collapsed) <- c(
    "gene_id", "seqnames", "start", "end", "strand", "gene_biotype", "gene_name"
  )
  collapsed$gene_name <- make.unique(names = collapsed$gene_name)
  gene.ranges <- makeGRangesFromDataFrame(
    df = collapsed,
    keep.extra.columns = TRUE
  )
  return(gene.ranges)
}


######### cicero functions ##########

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
