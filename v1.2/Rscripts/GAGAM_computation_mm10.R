########################################################################################
#
# INSTALL REQUIRED PACKAGES IF NEEDED
# 
# N.B. This may require the installation of local libraries. Please check the README
# file of the project for a list of required packages
#
########################################################################################
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
# GAGAM COMPUTATION FOR MOUSE MM10
# 
#######################################################################################
#
# NEEDED ELEMENTS FROM THE DATASETS PROCESSING SCRIPTS:
# 
# processed_ATAC_cds: THE SCATAC-SEQ DATA OBJECT
# labeled_peaks: LABELED PEAKS TABLE
# nmax
# brain: THE SEURAT OBJECT WITH THE ATAC DATA IN IT
# connection_table: THE CONNECTION OBJECT FROM THE CO-ACCESSIBILITY COMPUTATION
#
# OPTIONAL:
#
# CDS_RNA: THE OBJECT WITH PROCESSED SCRNA-SEQ, IF MULTI-OMIC DATASET
#
#######################################################################################
folder <- readline(prompt="Enter dataset folder name: ")
#paste0("../TMPResults/")"
#LOADING NECESSARY OBJECTS
processed_ATAC_cds <- readRDS(paste0("../TMPResults/Robjects/",folder,"/processed_ATAC_cds"))
labeled_peaks <- readRDS(paste0("../TMPResults/Robjects/",folder,"/labeled_peaks"))
nmax <- readRDS(paste0("../TMPResults/Robjects/",folder,"/nmax"))
brain <- readRDS(paste0("../TMPResults/Robjects/",folder,"/brain"))
connection_table <- readRDS(paste0("../TMPResults/Robjects/",folder,"/connection_table"))

# UNCOMMENT THE FOLLOWING LINE IF YOU WANT TO PERFORM COMPARATIVE ANALYSIS WITH SCRNA-SEQ DATA FROM A MULTI-OMIC DATASET
#CDS_RNA <- raedRDS(paste0("../TMPResults/Robjects/",folder,"/CDS_RNA"))

if (!(dir.exists("../TMPResults"))){
  dir.create("../TMPResults")
}

if (!(dir.exists("../TMPResults/classifications"))){
  dir.create("../TMPResults/classifications")
}
if (!(dir.exists(pasteo("../TMPResults/classifications/",folder)))){
  dir.create(pasteo("../TMPResults/classifications/",folder))
}


if (!(dir.exists("../TMPResults/GAM"))){
  dir.create("../TMPResults/GAM")
}
if (!(dir.exists(pasteo("../TMPResults/GAM/",folder)))){
  dir.create(pasteo("../TMPResults/GAM/",folder))
}


if (!(dir.exists("../TMPResults/metrics"))){
  dir.create("../TMPResults/metrics")
}
if (!(dir.exists(pasteo("../TMPResults/metrics/",folder)))){
  dir.create(pasteo("../TMPResults/metrics/",folder))
}

if (!(dir.exists("../TMPResults/IMAGES"))){
  dir.create("../TMPResults/IMAGES")
}
if (!(dir.exists(pasteo("../TMPResults/IMAGES/",folder)))){
  dir.create(pasteo("../TMPResults/IMAGES/",folder))
}

#####

#LOAD GENOME REFERENCE AND CHROMOSOMES INORMATION
refseq_anno <-  rtracklayer::readGFF("../DATA/Gene_2022/Genomes/mm10/refseq gene annotation/GCF_000001635.26_GRCm38.p6_genomic.gtf.gz")
chr2acc <- read.csv("../DATA/Gene_2022/Genomes/mm10/refseq gene annotation/chr2acc.txt", sep = "\t")

#SUBSET TO ONLY GENES
refseq_gene_anno <- refseq_anno[refseq_anno$seqid %in% chr2acc$Accession.version,]
refseq_gene_anno$seqid <- as.factor(as.character(refseq_gene_anno$seqid))
levels(refseq_gene_anno$seqid) <- chr2acc$X.Chromosome
refseq_gene_anno$seqid <- paste0("chr", refseq_gene_anno$seqid)
refseq_gene_anno <- refseq_gene_anno[refseq_gene_anno$type == "gene",]
refseq_gene_anno <- refseq_gene_anno[refseq_gene_anno$gene_biotype %in% c("protein_coding"),]

#DEFINING TSS AND FINDING TASS PEAKS
pos <- subset(refseq_gene_anno, strand == "+")
pos <- pos[order(pos$start),] 
pos$end <- pos$start + 1 

neg <- subset(refseq_gene_anno, strand == "-")
neg <- neg[order(neg$start, decreasing = TRUE),] 
neg$start <- neg$end - 1

refseq_gene_annotation_sub <- rbind(pos, neg)
refseq_gene_annotation_sub <- refseq_gene_annotation_sub[,c("seqid", "start", "end", "gene_id")]

names(refseq_gene_annotation_sub)[4] <- "gene"
processed_ATAC_cds <- annotate_cds_by_site(processed_ATAC_cds, refseq_gene_annotation_sub)


#CONSTRUCT PEAKS INFORMATION FROM LABELED PEAKS
peaks_multi_info <- fData(processed_ATAC_cds)
peaks_multi_info$site_name <- rownames(fData(processed_ATAC_cds))
peaks_multi_info$is.prom <- FALSE
peaks_multi_info$is.enhD <- FALSE
peaks_multi_info <- cbind(peaks_multi_info,labeled_peaks[7:(nmax+6)])

#PROM PEAKS
ppp <- peaks_multi_info
ppp <- as.data.frame(peaks_multi_info)
ppp[ppp == "prom"] <- TRUE

peaks_multi_info$is.prom <- apply(ppp, 1, any)
peaks_multi_info[is.na(peaks_multi_info$is.prom),]$is.prom <- FALSE
eaks_prom <- peaks_multi_info[peaks_multi_info$is.prom == "TRUE",]
prom_list <- rownames(peaks_prom)
non_prom_list <- rownames(peaks_multi_info)
non_prom_list <- setdiff(non_prom_list, prom_list)
non_prom_peaks <- peaks_multi_info[non_prom_list,]
non_prom_peaks[!is.na(non_prom_peaks$gene),]
peaks_prom[!is.na(peaks_prom$gene),]


#ENHD PEAKS
ppp <- peaks_multi_info
ppp <- as.data.frame(peaks_multi_info)
ppp[ppp == "enhD"] <- TRUE

peaks_multi_info$is.enhD <- FALSE
peaks_multi_info$is.enhD <- apply(ppp,1, any)
peaks_multi_info[is.na(peaks_multi_info$is.enhD),]$is.enhD <- FALSE

enhD_peaks_list <- rownames(peaks_multi_info[peaks_multi_info$is.enhD == TRUE,])
enhD_peaks_list <- setdiff(enhD_peaks_list, prom_list)


peaks_multi_info$PT <- !(is.na(peaks_multi_info$gene)) | as.vector(peaks_multi_info$is.prom)
####### PROMOTER CONTRIBUTION ########

refseq_peaks_prom <- peaks_multi_info[peaks_multi_info$PT == "TRUE",]
refseq_prom_list <- rownames(refseq_peaks_prom)

mm10 <- read.table("../DATA/Gene_2022/Genomes/mm10/mm10.chrom.sizes.txt")
mm10 <- mm10[1:21,]
mm10 <- Seqinfo(mm10[1:21,]$V1, seqlengths= mm10[1:21,]$V2)
mm10@genome[] <- "mm10"

refseq_GRanges <- makeGRangesFromDataFrame(refseq_gene_anno, seqinfo = mm10, seqnames.field = "seqid", keep.extra.columns = TRUE)

refseq_gene_closest_to_prom <- ClosestFeature(brain, gsub("_", "-", refseq_prom_list), annotation = refseq_GRanges)
refseq_gene_closest_to_prom <- refseq_gene_closest_to_prom[refseq_gene_closest_to_prom$distance <= 500,]

refseq_gene_closest_to_prom$query_region <- gsub("-","_",refseq_gene_closest_to_prom$query_region)

refseq_near_prom_list <- gsub("-","_",refseq_gene_closest_to_prom$query_region)


refseq_peaks_prom <-  refseq_peaks_prom[refseq_near_prom_list,]
refseq_peaks_prom$gene_anno <- refseq_gene_closest_to_prom$gene
refseq_peaks_prom$site_name <- rownames(refseq_peaks_prom)


refseq_promoter_peak_table <- refseq_peaks_prom[,c("gene","gene_anno","site_name")]
refseq_promoter_peak_table[is.na(refseq_promoter_peak_table$gene),]$gene <- refseq_promoter_peak_table[is.na(refseq_promoter_peak_table$gene),]$gene_anno
prom_gene_name <- levels(factor(refseq_promoter_peak_table$gene))

refseq_promoter_gene_mat <-
  Matrix::sparseMatrix(j=as.numeric(factor(refseq_promoter_peak_table$site_name)),
                       i=as.numeric(factor(refseq_promoter_peak_table$gene)),
                       x=1)

refseq_accessibility_mat <- exprs(processed_ATAC_cds)
refseq_accessibility_mat@x[refseq_accessibility_mat@x>0] <-1

refseq_promoter_activity_scores <- refseq_accessibility_mat[refseq_near_prom_list,, drop=FALSE]

colnames(refseq_promoter_gene_mat) = levels(factor(refseq_promoter_peak_table$site_name))
row.names(refseq_promoter_gene_mat) = levels(factor(refseq_promoter_peak_table$gene))
refseq_promoter_gene_mat <- refseq_promoter_gene_mat[,row.names(refseq_promoter_activity_scores)]

refseq_first_gene_matrix2 <- refseq_promoter_gene_mat %*% refseq_promoter_activity_scores
refseq_first_gene_matrix <- refseq_promoter_gene_mat %*% refseq_promoter_activity_scores
refseq_first_gene_matrix2 <- refseq_first_gene_matrix
refseq_first_gene_matrix2@x[refseq_first_gene_matrix2@x > 0] <- 1


####### EXON CONTRIBUTION #########

prom_list <- rownames(refseq_promoter_peak_table)
non_prom_list <- rownames(peaks_multi_info)
non_prom_list <- setdiff(non_prom_list, prom_list)
non_prom_list <- setdiff(non_prom_list, enhD_peaks_list)


refseq_exon_anno <- refseq_anno[refseq_anno$seqid %in% chr2acc$Accession.version,]
refseq_exon_anno$seqid <- as.factor(as.character(refseq_exon_anno$seqid))
levels(refseq_exon_anno$seqid) <- chr2acc$X.Chromosome
refseq_exon_anno$seqid <- paste0("chr", refseq_exon_anno$seqid)
refseq_exon_anno <- refseq_exon_anno[refseq_exon_anno$type == "exon",]
refseq_GRanges_exon <- makeGRangesFromDataFrame(refseq_exon_anno, seqinfo = mm10, seqnames.field = "seqid", keep.extra.columns = TRUE)

intragenetic_non_prom_peaks <- ClosestFeature(brain, gsub("_", "-", non_prom_list), refseq_GRanges_exon)
intragenetic_non_prom_peaks <- intragenetic_non_prom_peaks[intragenetic_non_prom_peaks$distance == 0,]
intragenetic_non_prom_peaks_list <- gsub("-","_",intragenetic_non_prom_peaks$query_region)

intragenetic_peaks <-  peaks_multi_info[intragenetic_non_prom_peaks_list,]
intragenetic_peaks$gene_anno <- intragenetic_non_prom_peaks$gene
intragenetic_peaks$site_name <- rownames(intragenetic_peaks)

intragenetic_peak_table <- intragenetic_peaks[,c("gene_anno","site_name")]

refseq_gene_annotation_sub <- refseq_gene_annotation_sub[refseq_gene_annotation_sub$gene %in% intragenetic_peak_table$gene_anno,]
refseq_gene_annotation_sub$TSS  <- paste0(refseq_gene_annotation_sub$seqid,"_",refseq_gene_annotation_sub$start,"_", refseq_gene_annotation_sub$end)

intragenetic_peak_table <- intragenetic_peak_table[intragenetic_peak_table$gene_anno %in% prom_gene_name,]

prova <- intragenetic_peak_table
prova <- as.data.frame(prova)
prova <- prova %>% rowwise() %>% mutate(TSS = refseq_gene_annotation_sub[refseq_gene_annotation_sub$gene == gene_anno,]$TSS) 
intragenetic_peak_table$TSS <- prova$TSS

if ("dist" %in% colnames(intragenetic_peak_table) == FALSE) {
  Peak1_cols <- split_peak_names(intragenetic_peak_table$site_name)
  Peak2_cols <- split_peak_names(intragenetic_peak_table$TSS)
  Peak1_bp <- round((as.integer(Peak1_cols[,3]) +
                       as.integer(Peak1_cols[,2])) / 2)
  Peak2_bp <- round((as.integer(Peak2_cols[,3]) +
                       as.integer(Peak2_cols[,2])) / 2)
  intragenetic_peak_table$dist <- abs(Peak2_bp - Peak1_bp)
}
intragenetic_peak_table$dist <- exp(-1*intragenetic_peak_table$dist/5000)


allgene_peaks_list <- rownames(intragenetic_peak_table)

refseq_allgene_gene_mat <-
  Matrix::sparseMatrix(j=as.numeric(factor(intragenetic_peak_table$site_name)),
                       i=as.numeric(factor(intragenetic_peak_table$gene_anno)),
                       x=intragenetic_peak_table$dist)

colnames(refseq_allgene_gene_mat) = levels(factor(intragenetic_peak_table$site_name))
row.names(refseq_allgene_gene_mat) = levels(factor(intragenetic_peak_table$gene_anno))
refseq_allgene_gene_mat <- refseq_allgene_gene_mat[row.names(refseq_allgene_gene_mat) %in% prom_gene_name,]



refseq_allgene_accessibility_mat <- exprs(processed_ATAC_cds)
refseq_allgene_accessibility_mat@x[refseq_allgene_accessibility_mat@x>0] <- 1
refseq_allgene_activity_scores <- refseq_allgene_accessibility_mat[allgene_peaks_list,, drop=FALSE]



refseq_allgene_gene_mat <- refseq_allgene_gene_mat[,row.names(refseq_allgene_activity_scores)]

refseq_allgene_gene_matrix <- refseq_allgene_gene_mat %*% refseq_allgene_activity_scores

refseq_allgene_gene_matrix <-  refseq_allgene_gene_matrix * refseq_first_gene_matrix2[rownames(refseq_allgene_gene_matrix),]


####### ENHD CONTRIBUTION ########
#### functions ####
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

#####
accessibility_mat <- exprs(processed_ATAC_cds)
accessibility_mat@x[accessibility_mat@x>0] <-1
#rownames(accessibility_mat) <- gsub(":", "_",rownames(accessibility_mat))
#rownames(accessibility_mat) <- gsub("-", "_",rownames(accessibility_mat))
con_val <- connection_table[connection_table$coaccess > 0,]
con_val <- con_val[!is.na(con_val$coaccess),]
coaccess <- signif(mean(con_val$coaccess), digits = 2)

if ("dist" %in% colnames(connection_table) == FALSE) {
  Peak1_cols <- split_peak_names(connection_table$Peak1)
  Peak2_cols <- split_peak_names(connection_table$Peak2)
  Peak1_bp <- round((as.integer(Peak1_cols[,3]) +
                       as.integer(Peak1_cols[,2])) / 2)
  Peak2_bp <- round((as.integer(Peak2_cols[,3]) +
                       as.integer(Peak2_cols[,2])) / 2)
  connection_table$dist <- abs(Peak2_bp - Peak1_bp)
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


promoter_peak_table <- refseq_promoter_peak_table[, c("site_name", "gene", "gene_anno")]



promoter_peak_table$site_name <- as.character(row.names(promoter_peak_table))
promoter_peak_table <- promoter_peak_table[!is.na(promoter_peak_table$gene),]
promoter_peak_table <- promoter_peak_table[,c("site_name", "gene")]
#promoter_peak_table$gene_anno <- as.character(promoter_peak_table$gene)

colnames(promoter_peak_table) <- c("peak", "gene")

dist_thresh=300000
coaccess
prom_enhD <- connection_table[(connection_table$Peak1 %in%
                      promoter_peak_table$peak &
                      connection_table$Peak2 %in%
                      enhD_peaks_list) | (connection_table$Peak2 %in%
                                            promoter_peak_table$peak &
                                            connection_table$Peak1 %in%
                                            enhD_peaks_list),]

#prom_enhD <- connection_table[(connection_table$Peak1 %in%
#                                 promoter_peak_table$peak &
#                                 connection_table$Peak2 %in%
#                                 enhD_peaks_list),]



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

promoter_conn_matrix <-
  prom_enhD_connectivity[unique(promoter_peak_table$peak),]

promoter_safe_sites <- intersect(rownames(promoter_conn_matrix),
                                 row.names(accessibility_mat))
distal_safe_sites <- intersect(colnames(promoter_conn_matrix),
                               row.names(accessibility_mat))
distal_safe_sites <- setdiff(distal_safe_sites, promoter_safe_sites)

promoter_access_mat_in_cicero_map <- accessibility_mat[promoter_safe_sites,, drop=FALSE]

distal_activity_scores <- accessibility_mat[distal_safe_sites,, drop=FALSE]

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
  promoter_conn_matrix[,distal_safe_sites, drop=FALSE] %*%
  scaled_site_weights
scaled_site_weights@x[scaled_site_weights@x > 1] <- 1

distal_activity_scores <- scaled_site_weights %*% distal_activity_scores

distal_activity_scores <- distal_activity_scores[row.names(promoter_access_mat_in_cicero_map),, drop=FALSE]

promoter_activity_scores <-  distal_activity_scores  + refseq_promoter_activity_scores 

promoter_gene_mat <-
  Matrix::sparseMatrix(j=as.numeric(factor(promoter_peak_table$peak)),
                       i=as.numeric(factor(promoter_peak_table$gene)),
                       x=1)
colnames(promoter_gene_mat) = levels(factor(promoter_peak_table$peak))
row.names(promoter_gene_mat) = levels(factor(promoter_peak_table$gene))
promoter_gene_mat <- promoter_gene_mat[,row.names(promoter_activity_scores)]
gene_activity_scores <- promoter_gene_mat %*% promoter_activity_scores

####### FINAL CONSTRUCTION ########

GAGAM2 <- gene_activity_scores 

# UNCOMMENT AND RUN THE FOLLOWING LINE ONLY IF ONE WANTS ALL THREE CONTRIBUTIONS
# GAGAM2[rownames(refseq_allgene_gene_matrix),] <- GAGAM2[rownames(refseq_allgene_gene_matrix),] + refseq_allgene_gene_matrix

GAGAM1 <- GAGAM2
final_GAM <- GAGAM1

num_genes <- pData(processed_ATAC_cds)$num_genes_expressed
names(num_genes) <- row.names(pData(processed_ATAC_cds))

# NORMALIZATION OF GAGAM
final_GAM <- normalize_gene_activities(prova, num_genes)


####### GAGAM ANALYSIS #######
first_cell <- colnames(final_GAM)
lenght1 <- length(first_cell)
first_gene <- row.names(final_GAM)
lenght2 <- length(first_gene)
c_c <- matrix(first_cell, nrow = lenght1, dimnames = list(first_cell,c("Cells")))
c_g<- matrix(first_gene, nrow = lenght2, dimnames = list(first_gene,c("gene_short_name")))


final_GAM <- final_GAM[!Matrix::rowSums(final_GAM) == 0, 
                               !Matrix::colSums(final_GAM) == 0]
final_GAM_cds <-  suppressWarnings(new_cell_data_set(final_GAM, cell_metadata = c_c, gene_metadata = c_g))
rowData(final_GAM_cds)$gene_short_name <- rownames(final_GAM_cds)
final_GAM_cds <- detect_genes(final_GAM_cds)
final_GAM_cds <- estimate_size_factors(final_GAM_cds)
final_GAM_cds <- preprocess_cds(final_GAM_cds, method = "PCA",num_dim = 50)

final_GAM_cds <- reduce_dimension(final_GAM_cds, reduction_method = 'UMAP', 
                                  preprocess_method = "PCA")

#SET YOUR OWN RESOLUTION
final_GAM_cds = cluster_cells(final_GAM_cds, resolution=1.2e-3)

pdf(pasteo("../TMPResults/IMAGES/",folder,"/GAGAM_plot.pdf"))
plot_cells(final_GAM_cds, cell_size = 1.2, group_label_size = 5)
dev.off()

#PUT INSIDE THE cds THE ATAC CLUSTERING
class <- as.data.frame(processed_ATAC_cds@clusters@listData[["UMAP"]][["clusters"]])
colnames(class) <- "CLASS"
final_GAM_cds@colData@listData[["label"]] <- class$CLASS

#plot_cells(final_GAM_cds, color_cells_by = "label", label_groups_by_cluster = FALSE ,cell_size = 1.2, group_label_size = 5, label_cell_groups = FALSE)

#CLUSTERING RESULTS
final_GAM_first <- as.data.frame(final_GAM_cds@clusters@listData[["UMAP"]][["clusters"]])
colnames(final_GAM_first) <- "CLASS"


#METRICS
ARI <-ARI(class$CLASS, final_GAM_first$CLASS)
AMI <- AMI(class$CLASS, final_GAM_first$CLASS)
ARI
AMI

#housekeeping_genes <- read.table("DATA/buenrostro/HK_genes.txt")
#housekeeping_genes <- housekeeping_genes$V1
#diff <- setdiff(housekeeping_genes, prom_gene_name)
#our_house <- setdiff(housekeeping_genes, diff)
#write.table(our_house, file='housekeeping_genes.txt', quote=FALSE, sep='\t')

#SAVING GAGAM AND CLUSTERING
tsv <- data.matrix(exprs(final_GAM_cds))

write.table(tsv, file=paste0('../TMPResults/GAM/',folder,'/gagam.tsv'), quote=FALSE, sep='\t', col.names = NA)
write.table(final_GAM_first, file=paste0('../TMPREsults/classifications/',folder,'/gagam_classification.tsv'), quote=FALSE, sep='\t', col.names = NA)
#write.table(class, file=paste0('../TMPResults/classifications/',folder,'atac_classification.tsv'), quote=FALSE, sep='\t', col.names = NA)


####### SCRNA-SEQ CONFRONTATION #######

if (exists("CDS_RNA")) {
final_GAM_cds@colData@listData[["cluster"]] <- final_GAM_first$CLASS

markers <- top_markers(final_GAM_cds, group_cells_by = "label")
top_specific_markers_gagam <- markers %>%
  filter(specificity >= 0.15) %>%
  group_by(cell_group) %>%
  top_n(1, marker_score)
top_specific_marker_gagam_ids <- unique(top_specific_markers_gagam %>% pull(gene_id))
# uncommnet to plot the DE features
#plot_genes_by_group(final_GAM_cds,
#                    top_specific_marker_gagam_ids,
#                    group_cells_by="cluster",
#                    ordering_type="maximal_on_diag")


S <- as.Seurat(cds, counts = "counts", data = NULL)
CDS_RNA@colData@listData[["activity"]] <- first$CLASS
C <- as.Seurat(CDS_RNA, counts = "counts", data = NULL)
C<- NormalizeData(C)
C <- ScaleData(C)
S<- NormalizeData(S)
S <- ScaleData(S)

pdf(paste0("../TMPResults/IMAGES/",folder,"/vlnplot_activity.pdf"))
VlnPlot(S, features = c("MS4A1", "AGAP1", "KLF4"), group.by = "cluster")
dev.off()
pdf(paste0("../TMPResults/IMAGES/",folder,"vlnplot_expression.pdf"))
VlnPlot(C, features = c("MS4A1", "AGAP1", "KLF4"), group.by = "activity")
dev.off()

Idents(object = S) <- "cluster"
Idents(object = C) <- "activity"

pdf(paste0("../TMPResults/IMAGES/",folder,"heatmap_activity.pdf"))
DoHeatmap(S, features = top_specific_marker_gagam_ids, size = 3)
dev.off()
pdf(paste0("../TMPResults/IMAGES/",folder,"heatmap_expression.pdf")
DoHeatmap(C, features = top_specific_marker_gagam_ids, size = 3)
dev.off()

}

#######################################################################################
#
# END
# 
#######################################################################################



