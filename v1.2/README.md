# GAGAM v1.2

GAGAM is Genomic Annotated Gene Activity Matrix, a model-driven approach to interpret scATAC-seq accessibility data, to obtain information about gene activity. It is a fundamental step to correlate scATAC-seq with ScRNA-seq. In other words, to link accessibility and expression of genes.

GAGAM v1.0 was used in the experimental setup published in:
Martini, L., Bardini, R., Savino, A., Di Carlo, S. (2022). GAGAM: A Genomic Annotation-Based Enrichment of scATAC-seq Data for Gene Activity Matrix. In: Rojas, I., Valenzuela, O., Rojas, F., Herrera, L.J., Ortuño, F. (eds) Bioinformatics and Biomedical Engineering. IWBBIO 2022. Lecture Notes in Computer Science(), vol 13347. Springer, Cham. https://doi.org/10.1007/978-3-031-07802-6_2.

GAGAM v1.2 was used in the experimental setup in:
Martini, L., Bardini, R., Savino, A., Di Carlo, S. (2022). GGAGAM v1.2: an improvement on peak labeling and Genomic Annotated Gene Activity Matrix construction. Submmitted for revision for MDPI Genes journal.

This README file lists all available scripts and describes the full workflow to reproduce the experiments reported in the paper. 

They are divided into:

1. R scripts
1. Python scripts

## R scripts

1. **10x_v1_PBMChuman_processing.R**: contains the scripts to process the raw scATAC-seq data,to compute Cicero and GeneScoring, and process them, and prepare the objects need for GAGAM computation, for the 10X V1.0.1 PBMC dataset, containing human PBMCs.

2. **10X_V2_PBMC_processing.R**: contains the scripts to process the raw scATAC-seq data,to compute Cicero and GeneScoring, and process them, and prepare the objects need for GAGAM computation, for the 10X V2.0.0 PBMC dataset, containing human PBMCs.

3. **buenrostro_processing.R** : ccontains the scripts to process the raw scATAC-seq data,to compute Cicero and GeneScoring, and process them, and prepare the objects need for GAGAM computation, for the Buenrostro2018 dataset, containing bone marrow cells.

4. **10x_multiome_processing.R**: contains the scripts to process the raw scATAC-seq data,to compute Cicero and GeneScoring, and process them, and prepare the objects need for GAGAM computation, for the 10X Multiome Controller PBMC dataset, containing mouse brain cells.

5. **10x_multiome_X_processing.R**: contains the scripts to process the raw scATAC-seq data,to compute Cicero and GeneScoring, and process them, and prepare the objects need for GAGAM computation, for the 10X Multiome Chromium X PBMC dataset, containing mouse brain cells.

6. **10x_v2_mousebrain_processing.R**: contains the scripts to process the raw scATAC-seq data,to compute Cicero and GeneScoring, and process them, and prepare the objects need for GAGAM computation, for the 10X V2.1.0 Brain dataset, containing mouse brain cells.

7. **GAGAM_compuation_hg38**: contains the scripts for the GAGAM computation and its processing, in case of data alligned to the human hg38 genome assembly. Optionally it perform comparative analysis with scRNA-seq data, in case of multi-omic datasets.

8. **GAGAM_computation_mm10**: contains the scripts for the GAGAM computation and its processing, in case of data alligned to the mouse mm10 genome assembly. Optionally it perform comparative analysis with scRNA-seq data, in case of multi-omic datasets.

## Python scripts

Python scripts use the python3 language, and include:

1. **classifyPeaks.py**: labels peaks from a peak file using cCRECombined genomic annotation tracks;
2. **computeRAGI.py**: computes the RAGI scores for a GAM versus another classification, leveraging a list of housekeeping genes and a list of marker genes.

# Procedure

## Phase 0: Data preparation

For each dataset, download the data from the sources

- 10X v1.0.1 PBMC: https://www.10xgenomics.com/resources/datasets/5-k-peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-1-standard-1-0-1
- 10X V2.0.0 PBMC: https://www.10xgenomics.com/resources/datasets/5-k-peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-next-gem-v-1-1-1-1-standard-2-0-0
- Buenrostro2018: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96769
- 10X Multiome Controller PBMC: https://www.10xgenomics.com/resources/datasets/10-k-human-pbm-cs-multiome-v-1-0-chromium-controller-1-standard-2-0-0
- 10X Multiome Chromium X PBMC: https://www.10xgenomics.com/resources/datasets/10-k-human-pbm-cs-multiome-v-1-0-chromium-x-1-standard-2-0-0
- 10X V2.1.0 Brain: https://www.10xgenomics.com/resources/datasets/8k-adult-mouse-cortex-cells-atac-v1-1-chromium-x-1-1-standard

Organize the data following the proposed folder structure in the TMPDATA folder:

- 10x_V1_PBMChuman
  - barcodes.tsv
  - matrix.mtx
  - peaks.bed
  - atac_v1_pbmc_5k_filtered_peak_bc_matrix.h5
- 10X_V2_PBMC
  - barcodes.tsv
  - matrix.mtx
  - peaks.bed
  - atac_pbmc_5k_nextgem_filtered_peak_bc_matrix.h5
- buenrostro
  - GSE96769_PeakFile_20160207.csv.bed
  - GSE96769_scATACseq_counts.txt
  - hglft_genome_241cf_9a53c0.bed
  - metadata.tsv
- 10x_PBMC_Multiome_Controller
  - barcodes.tsv
  - matrix.mtx
  - features.bed
- 10x_PBMC_Multiome_ChromiumX
  - barcodes.tsv
  - matrix.mtx
  - features.bed
- 10x_v2_mousebrain
  - barcodes.tsv
  - matrix.mtx
  - peaks.bed

## Phase 1: Peak labeling

The first step is to label the peaks with Genomic annotation tracks.

## Input

- File containig peaks for each dataset
- Species appropriate cCRE Combined annotation tracks

## Script

- **classifyPeaks.py**

## Syntax

Execute the script from the v1.2 folder with the syntax:
```
##10X V1.0.1 PBMC peak labeling
python3 classifyPeaks.py '../TMPDATA/10x_V1_PBMChuman/peaks.bed' '../DATA/Genomic_annotations/encodeCcreCombined_human.bb'

##10X V2.0.0 PBMC peak labeling
python3 classifyPeaks.py '../TMPDATA/10X_V2_PBMC/peaks.bed' '../DATA/Genomic_annotations/encodeCcreCombined_human.bb'

##Buenrostro2018 peak labeling
python3 classifyPeaks.py '../TMPDATA/buenrostro/GSE96769_PeakFile_20160207.csv.bed' '../DATA/Genomic_annotations/encodeCcreCombined_human.bb'

##10X Multiome Controller PBMC peak labeling
python3 classifyPeaks.py '../TMPDATA/10x_PBMC_Multiome_Controller/features.bed' '../DATA/Genomic_annotations/encodeCcreCombined_human.bb'

##10X Multiome Chromium X PBMC peak labeling
python3 classifyPeaks.py '../TMPDATA/10x_PBMC_Multiome_ChromiumX/features.bed' '../DATA/Genomic_annotations/encodeCcreCombined_human.bb'

##10X V2.1.0 Brain peak labeling
python3 classifyPeaks.py '../TMPDATA/10x_v2_mousebrain/peaks.bed'../DATA/Genomic_annotations/encodeCcreCombined_mouse.bb'

```

## Output

The script will output a file containing the labeled peaks for each dataset in the respective folder:

- ../TMPResults/labeled_peaks/10x_V1_PBMChuman/cCRE_labeled_peaks.tsv
- ../TMPResults/labeled_peaks/10X_V2_PBMC/cCRE_labeled_peaks.tsv
- ../TMPResults/labeled_peaks/buenrostro/cCRE_labeled_peaks.tsv
- ../TMPResults/labeled_peaks/10x_PBMC_Multiome_Controller/cCRE_labeled_peaks.tsv
- ../TMPResults/labeled_peaks/10x_PBMC_Multiome_ChromiumX/cCRE_labeled_peaks.tsv
- ../TMPResults/labeled_peaks/10x_v2_mousebrain/cCRE_labeled_peaks.tsv

## Phase 2: Dataset processing

The second phase consists in scATAC-seq data standard processing, construction of Cicero and GeneScoring GAM, and clustering classification of all of them. Optionally, it takes and process also the scRNA-seq data from the multi-omic datasets, and create the objects for the comparaive analysis with gene expression.

### Input

Path for 10X V1.0.1 PBMC inputs:
  - ../TMPResults/leabeled_peaks/10x_V1_PBMChuman/cCRE_labeled_peaks.csv
  - ../TMPDATA/10x_V1_PBMChuman/
  - ../DATA/Genomes/hg38/GCF_000001405.39_GRCh38.p13_genomic.gtf

Path for 10X V2.0.0 PBMC inputs:
  - ../TMPResults/labeled_peaks/10X_V2_PBMC/cCRE_labeled_peaks.csv
  - ../TMPDATA/10X_V2_PBMC/
  - ../DATA/Genomes/hg38/GCF_000001405.39_GRCh38.p13_genomic.gtf

Path for Buenrostro2018 inputs:
  - ../TMPResults/labeled_peaks/buenrostro/cCRE_labeled_peaks.csv
  - ../TMPDATA/buenrostro/
  - ../DATA/Genomes/hg38/GCF_000001405.39_GRCh38.p13_genomic.gtf

Path for 10X Multiome Controller PBMC inputs:
  - ../TMPResults/labeled_peaks/
  - /cCRE_labeled_peaks.csv
  - ../TMPDATA/10x_PBMC_Multiome_Controller/
  - ../DATA/Genomes/hg38/GCF_000001405.39_GRCh38.p13_genomic.gtf

Path for 10X Multiome Chromium X PBMC inputs:
  - ../TMPResults/labeled_peaks/10x_PBMC_Multiome_ChromiumX/cCRE_labeled_peaks.csv
  - ../TMPDATA/10x_PBMC_Multiome_ChromiumX/
  - ../DATA/Genomes/mm10/GCF_000001635.26_GRCm38.p6_genomic.gtf

Path for 10X V2.1.0 Brain inputs:
  - ../TMPResults/labeled_peaks/10x_v2_mousebrain/cCRE_labeled_peaks.csv
  - ../TMPDATA/10x_v2_mousebrain/
  - ../DATA/Genomes/mm10/GCF_000001635.26_GRCm38.p6_genomic.gtf


### Script

- **10x_v1_PBMChuman_processing.R** 
- **10X_V2_PBMC_processing.R**
- **buenrostro_processing.R**
- **10x_multiome_processing.R**
- **10x_multiome_X_processing.R**
- **10x_v2_mousebrain_processing.R**


### Syntax

Execute the script from the v1.2 folder with the syntax:
```
##10X V1.0.1 PBMC processing
Rscript 10x_v1_PBMChuman_processing.R

##10X V2.0.0 PBMC processing
Rscript 10X_V2_PBMC_processing.R

##Buenrostro2018 processing
Rscript buenrostro_processing.R 

##10X Multiome Controller PBMC processing
Rscript 10x_multiome_processing.R 

##10X Multiome Chromium X PBMC processing
Rscript 10x_multiome_X_processing.R

##10X V2.1.0 Brain processing
Rscript 10x_v2_mousebrain_processing.R

```

### Output

Path for 10x_v1_PBMChuman_processing.R outputs:
  - ../TMPResults/GAM/10x_V1_PBMChuman/cicero.tsv
  - ../TMPResults/GAM/10x_V1_PBMChuman/genescoring.tsv
  - ../TMPResults/classifications/10x_V1_PBMChuman/cicero_classifications.tsv
  - ../TMPResults/classifications/10x_V1_PBMChuman/genescoring_classifications.tsv
  - ../TMPResults/classifications/10x_V1_PBMChuman/atac_classification.tsv
  - ../TMPResults/Robjects/10x_V1_PBMChuman/processed_ATAC_cds
  - ../TMPResults/Robjects/10x_V1_PBMChuman/connection_table

Path for 10X_V2_PBMC_processing.R outputs:
  - ../TMPResults/GAM/10X_V2_PBMC/cicero.tsv
  - ../TMPResults/GAM/10X_V2_PBMC/genescoring.tsv
  - ../TMPResults/classifications/10X_V2_PBMC/cicero_classifications.tsv
  - ../TMPResults/classifications/10X_V2_PBMC/genescoring_classifications.tsv
  - ../TMPResults/classifications/10X_V2_PBMC/atac_classification.tsv
  - ../TMPResults/Robjects/10X_V2_PBMC/processed_ATAC_cds
  - ../TMPResults/Robjects/10X_V2_PBMC/connection_table

Path for buenrostro_processing.R outputs:
  - ../TMPResults/GAM/buenrostro/cicero.tsv
  - ../TMPResults/GAM/buenrostro/genescoring.tsv
  - ../TMPResults/classifications/buenrostro/cicero_classifications.tsv
  - ../TMPResults/classifications/buenrostro/genescoring_classifications.tsv
  - ../TMPResults/classifications/buenrostro/atac_classification.tsv
  - ../TMPResults/Robjects/buenrostro/processed_ATAC_cds
  - ../TMPResults/Robjects/buenrostro/connection_table

Path for 10x_multiome_processing.R outputs:
  - ../TMPResults/GAM/10x_PBMC_Multiome_Controller/cicero.tsv
  - ../TMPResults/GAM/10x_PBMC_Multiome_Controller/genescoring.tsv
  - ../TMPResults/classifications/10x_PBMC_Multiome_Controller/cicero_classifications.tsv
  - ../TMPResults/classifications/10x_PBMC_Multiome_Controller/genescoring_classifications.tsv
  - ../TMPResults/classifications/10x_PBMC_Multiome_Controller/atac_classification.tsv
  - ../TMPResults/Robjects/10x_PBMC_Multiome_Controller/processed_ATAC_cds
  - ../TMPResults/Robjects/10x_PBMC_Multiome_Controller/connection_table

Path for 10x_multiome_X_processing.R outputs:
  - ../TMPResults/GAM/10x_PBMC_Multiome_ChromiumX/cicero.tsv
  - ../TMPResults/GAM/10x_PBMC_Multiome_ChromiumX/genescoring.tsv
  - ../TMPResults/classifications/10x_PBMC_Multiome_ChromiumX/cicero_classifications.tsv
  - ../TMPResults/classifications/10x_PBMC_Multiome_ChromiumX/genescoring_classifications.tsv
  - ../TMPResults/classifications/10x_PBMC_Multiome_ChromiumX/atac_classification.tsv
  - ../TMPResults/Robjects/10x_PBMC_Multiome_ChromiumX/processed_ATAC_cds
  - ../TMPResults/Robjects/10x_PBMC_Multiome_ChromiumX/connection_table

Path for 10x_v2_mousebrain_processing.R outputs:
  - ../TMPResults/GAM/10x_v2_mousebrain/cicero.tsv
  - ../TMPResults/GAM/10x_v2_mousebrain/genescoring.tsv
  - ../TMPResults/classifications/10x_v2_mousebrain/cicero_classifications.tsv
  - ../TMPResults/classifications/10x_v2_mousebrain/genescoring_classifications.tsv
  - ../TMPResults/classifications/10x_v2_mousebrain/atac_classification.tsv
  - ../TMPResults/Robjects/10x_v2_mousebrain/processed_ATAC_cds
  - ../TMPResults/Robjects/10x_v2_mousebrain/connection_table



## Phase 3: GAGAM construction

### Input

The inputs are all the Robjects saved by the processing scripts in the respective folders.Moreover, at the start of the execution it will ask to provide the dataset folder name. This will allow the script to find the rigth inputs and save the outputs in the respective folders. Optionally, perform the comparative analysis between gene activity and gene expression, for the multi-omic datasets.

Path for 10X V1.0.1 PBMC inputs:
  - ../TMPResults/Robjects/10x_V1_PBMChuman/

Path for 10X V2.0.0 PBMC inputs:
  - ../TMPResults/Robjects/10X_V2_PBMC/

Path for Buenrostro2018 inputs:
  - ../TMPResults/Robjects/buenrostro/

Path for 10X Multiome Controller PBMC inputs:
  - ../TMPResults/Robjects/10x_PBMC_Multiome_Controller/

Path for 10X Multiome Chromium X PBMC inputs:
  - ../TMPResults/Robjects/10x_PBMC_Multiome_ChromiumX/

Path for 10X V2.1.0 Brain inputs:
  - ../TMPResults/Robjects/10x_v2_mousebrain/

### Script

- **GAGAM_computation_hg38**
- **GAGAM_computation_mm10**

### Syntax

Execute the script from the v1.2 folder with the syntax:

```
##10X V1.0.1 PBMC GAGAM construction
Rscript GAGAM_computation_hg38.R

##10X V2.0.0 PBMC GAGAM construction
Rscript GAGAM_computation_hg38.R

##Buenrostro2018 GAGAM construction
Rscript GAGAM_computation_hg38.R 

##10X Multiome Controller PBMC GAGAM construction
Rscript GAGAM_computation_hg38.R 

##10X Multiome Chromium X PBMC GAGAM construction
Rscript GAGAM_computation_hg38.R

##10X V2.1.0 Brain GAGAM construction
Rscript GAGAM_computation_mm10.R
```

### Output

The script execution outputs the GAGAM, its clustering classifications, IMAGES, and ARI and AMI tables. 

Path for 10X V1.0.1 PBMC outputs:
  - ../TMPResults/GAM/10x_V1_PBMChuman/gagam.tsv
  - ../TMPResults/classifications/10x_V1_PBMChuman/gagam_classifications.tsv
  - ../TMPResults/metrics/10x_V1_PBMChuman/ARI_AMI_table.tsv
  - ../TMPResults/IMAGES/10x_V1_PBMChuman/

Path for 10X V2.0.0 PBMC outputs:
  - ../TMPResults/GAM/10X_V2_PBMC/gagam.tsv
  - ../TMPResults/classifications/10X_V2_PBMC/gagam_classifications.tsv
  - ../TMPResults/metrics/10X_V2_PBMC/ARI_AMI_table.tsv
  - ../TMPResults/IMAGES/10X_V2_PBMC/

Path for Buenrostro2018 outputs:
  - ../TMPResults/GAM/buenrostro/gagam.tsv
  - ../TMPResults/classifications/buenrostro/gagam_classifications.tsv
  - ../TMPResults/metrics/buenrostro/ARI_AMI_table.tsv
  - ../TMPResults/IMAGES/buenrostro/

Path for 10X Multiome Controller PBMC outputs:
  - ../TMPResults/GAM/10x_PBMC_Multiome_Controller/gagam.tsv
  - ../TMPResults/classifications/10x_PBMC_Multiome_Controller/gagam_classifications.tsv
  - ../TMPResults/metrics/10x_PBMC_Multiome_Controller/ARI_AMI_table.tsv
  - ../TMPResults/IMAGES/10x_PBMC_Multiome_Controller/

Path for 10X Multiome Chromium X PBMC outputs:
  - ../TMPResults/GAM/10x_PBMC_Multiome_ChromiumX/gagam.tsv
  - ../TMPResults/classifications/10x_PBMC_Multiome_ChromiumX/gagam_classifications.tsv
  - ../TMPResults/metrics/10x_PBMC_Multiome_ChromiumX/ARI_AMI_table.tsv
  - ../TMPResults/IMAGES/10x_PBMC_Multiome_ChromiumX/

Path for 10X V2.1.0 Brain outputs:
  - ../TMPResults/GAM/10x_v2_mousebrain/gagam.tsv
  - ../TMPResults/classifications/10x_v2_mousebrain/gagam_classifications.tsv
  - ../TMPResults/metrics/10x_v2_mousebrain/ARI_AMI_table.tsv
  - ../TMPResults/IMAGES/10x_v2_mousebrain/

## Phase 4: Metrics Calculation

The forth phase is the calculation of RAGI scores between each GAM and all the classifications for each dataset. 

### Input

For each dataset, the script takes as inputs the file for the GAM under analysis, the path of the folder containing all classification results to perform the comparison, and paths to species- and tissue-specific housekeeping and marker genes list files.

### Script

- **computeRAGI.py**

### Syntax

Compute RAGI for each GAM over classification results, for each dataset. Execute the script from the v1.2 folder with the syntax:

```
#### 10X V1.0.1 PBMC ######
#GAGAM v1.2
python3 computeRAGI.py '../TMPResults/GAM/10x_V1_PBMChuman/gagam.tsv' '../TMPResults/classifications/10x_V1_PBMChuman/' '../DATA/Genes_2022/Housekeeping_genes/Housekeeping_genes_human.txt' '../DATA/Genes_2022/Marker_genes/human_PBMC_gene_markers.csv'
#Cicero GAM
python3 computeRAGI.py '../TMPResults/GAM/10x_V1_PBMChuman/cicero.tsv' '../TMPResults/classifications/10x_V1_PBMChuman/' '../DATA/Genes_2022/Housekeeping_genes/Housekeeping_genes_human.txt' '../DATA/Genes_2022/Marker_genes/human_PBMC_gene_markers.csv'
#GeneScoring GAM
python3 computeRAGI.py '../TMPResults/GAM/10x_V1_PBMChuman/genescoring.tsv' '../TMPResults/classifications/10x_V1_PBMChuman/' '../DATA/Genes_2022/Housekeeping_genes/Housekeeping_genes_human.txt' '../DATA/Genes_2022/Marker_genes/human_PBMC_gene_markers.csv'
```
```
#### 10X V2.0.0 PBMC ######
#GAGAM v1.2
python3 computeRAGI.py '../TMPResults/GAM/10X_V2_PBMC/gagam.tsv' '../TMPResults/classifications/10X_V2_PBMC/' '../DATA/Genes_2022/Housekeeping_genes/Housekeeping_genes_human.txt' '../DATA/Genes_2022/Marker_genes/human_PBMC_gene_markers.csv'
#Cicero GAM
python3 computeRAGI.py '../TMPResults/GAM/10X_V2_PBMC/cicero.tsv' '../TMPResults/classifications/10X_V2_PBMC/' '../DATA/Genes_2022/Housekeeping_genes/Housekeeping_genes_human.txt' '../DATA/Genes_2022/Marker_genes/human_PBMC_gene_markers.csv'
#GeneScoring GAM
python3 computeRAGI.py '../TMPResults/GAM/10X_V2_PBMC/genescoring.tsv' '../TMPResults/classifications/10X_V2_PBMC/' '../DATA/Genes_2022/Housekeeping_genes/Housekeeping_genes_human.txt' '../DATA/Genes_2022/Marker_genes/human_PBMC_gene_markers.csv'
```
```
#### Buenrostro2018 ######
#GAGAM v1.2
python3 computeRAGI.py '../TMPResults/GAM/buenrostro/gagam.tsv' '../TMPResults/classifications/buenrostro/' '../DATA/Genes_2022/Housekeeping_genes/Housekeeping_genes_human.txt' '../DATA/Genes_2022/Marker_genes/human_PBMC_and_bone_marrow_gene_markers.csv'
#Cicero GAM
python3 computeRAGI.py '../TMPResults/GAM/buenrostro/cicero.tsv' '../TMPResults/classifications/buenrostro/' '../DATA/Genes_2022/Housekeeping_genes/Housekeeping_genes_human.txt' '../DATA/Genes_2022/Marker_genes/human_PBMC_and_bone_marrow_gene_markers.csv'
#GeneScoring GAM
python3 computeRAGI.py '../TMPResults/GAM/buenrostro/genescoring.tsv' '../TMPResults/classifications/buenrostro/' '../DATA/Genes_2022/Housekeeping_genes/Housekeeping_genes_human.txt' '../DATA/Genes_2022/Marker_genes/human_PBMC_and_bone_marrow_gene_markers.csv'
```
```
#### 10X Multiome Controller PBMC ######
#GAGAM v1.2
python3 computeRAGI.py '../TMPResults/GAM/10x_PBMC_Multiome_Controller/gagam.tsv' '../TMPResults/classifications/10x_PBMC_Multiome_Controller/' '../DATA/Genes_2022/Housekeeping_genes/Housekeeping_genes_human.txt' '../DATA/Genes_2022/Marker_genes/human_PBMC_Pliner2019_gene_markers.csv'
#Cicero GAM
python3 computeRAGI.py '../TMPResults/GAM/10x_PBMC_Multiome_Controller/cicero.tsv' '../TMPResults/classifications/10x_PBMC_Multiome_Controller/' '../DATA/Genes_2022/Housekeeping_genes/Housekeeping_genes_human.txt' '../DATA/Genes_2022/Marker_genes/human_PBMC_Pliner2019_gene_markers.csv'
#GeneScoring GAM
python3 computeRAGI.py '../TMPResults/GAM/10x_PBMC_Multiome_Controller/genescoring.tsv' '../TMPResults/classifications/10x_PBMC_Multiome_Controller/' '../DATA/Genes_2022/Housekeeping_genes/Housekeeping_genes_human.txt' '../DATA/Genes_2022/Marker_genes/human_PBMC_Pliner2019_gene_markers.csv'
```
```
#### 10X Multiome Chromium X PBMC ######
#GAGAM v1.2
python3 computeRAGI.py '../TMPResults/GAM/10x_PBMC_Multiome_ChromiumX/gagam.tsv' '../TMPResults/classifications/10x_PBMC_Multiome_ChromiumX/' '../DATA/Genes_2022/Housekeeping_genes/Housekeeping_genes_human.txt' '../DATA/Genes_2022/Marker_genes/human_PBMC_Pliner2019_gene_markers.csv'
#Cicero GAM
python3 computeRAGI.py '../TMPResults/GAM/10x_PBMC_Multiome_ChromiumX/cicero.tsv' '../TMPResults/classifications/10x_PBMC_Multiome_ChromiumX/' '../DATA/Genes_2022/Housekeeping_genes/Housekeeping_genes_human.txt' '../DATA/Genes_2022/Marker_genes/human_PBMC_Pliner2019_gene_markers.csv'
#GeneScoring GAM
python3 computeRAGI.py '../TMPResults/GAM/10x_PBMC_Multiome_ChromiumX/genescoring.tsv' '../TMPResults/classifications/10x_PBMC_Multiome_ChromiumX/' '../DATA/Genes_2022/Housekeeping_genes/Housekeeping_genes_human.txt' '../DATA/Genes_2022/Marker_genes/human_PBMC_Pliner2019_gene_markers.csv'
```

```
#### 10X V2.1.0 Brain ######
#GAGAM v1.2
python3 computeRAGI.py '../TMPResults/GAM/10x_v2_mousebrain/gagam.tsv' '../TMPResults/classifications/10x_v2_mousebrain/' '../DATA/Genes_2022/Housekeeping_genes/Housekeeping_genes_mouse.txt' '../DATA/Genes_2022/Marker_genes/mouse_brain_gene_markers.csv'
#Cicero GAM
python3 computeRAGI.py '../TMPResults/GAM/10x_v2_mousebrain/cicero.tsv' '../TMPResults/classifications/10x_v2_mousebrain/' '../DATA/Genes_2022/Housekeeping_genes/Housekeeping_genes_mouse.txt' '../DATA/Genes_2022/Marker_genes/mouse_brain_gene_markers.csv'
#GeneScoring GAM
python3 computeRAGI.py '../TMPResults/GAM/10x_v2_mousebrain/genescoring.tsv' '../TMPResults/classifications/10x_v2_mousebrain/' '../DATA/Genes_2022/Housekeeping_genes/Housekeeping_genes_mouse.txt' '../DATA/Genes_2022/Marker_genes/mouse_brain_gene_markers.csv'
```

### Output

Path for 10x_v1_PBMChuman RAGI outputs:
  - ../TMPResults/metrics/10x_v1_PBMChuman/gagam_RAGI_table.csv
  - ../TMPResults/metrics/10x_v1_PBMChuman/cicero_RAGI_table.csv
  - ../TMPResults/metrics/10x_v1_PBMChuman/genescoring_RAGI_table.csv

Path for 10X_V2_PBMC outputs:
  - ../TMPResults/metrics/10X_V2_PBMC/gagam_RAGI_table.csv
  - ../TMPResults/metrics/10X_V2_PBMC/cicero_RAGI_table.csv
  - ../TMPResults/metrics/10X_V2_PBMC/genescoring_RAGI_table.csv

Path for buenrostro  outputs:
  - ../TMPResults/metrics/buenrostro/gagam_RAGI_table.csv
  - ../TMPResults/metrics/buenrostro/cicero_RAGI_table.csv
  - ../TMPResults/metrics/buenrostro/genescoring_RAGI_table.csv

Path for 10x_PBMC_Multiome_Controller  outputs:
  - ../TMPResults/metrics/10x_PBMC_Multiome_Controller/gagam_RAGI_table.csv
  - ../TMPResults/metrics/10x_PBMC_Multiome_Controller/cicero_RAGI_table.csv
  - ../TMPResults/metrics/10x_PBMC_Multiome_Controller/genescoring_RAGI_table.csv

Path for 10x_PBMC_Multiome_ChromiumX outputs:
  - ../TMPResults/metrics/10x_PBMC_Multiome_ChromiumX/gagam_RAGI_table.csv
  - ../TMPResults/metrics/10x_PBMC_Multiome_ChromiumX/cicero_RAGI_table.csv
  - ../TMPResults/metrics/10x_PBMC_Multiome_ChromiumX/genescoring_RAGI_table.csv

Path for 10x_v2_mousebrain outputs:
  - ../TMPResults/metrics/10x_v2_mousebrain/gagam_RAGI_table.csv
  - ../TMPResults/metrics/10x_v2_mousebrain/cicero_RAGI_table.csv
  - ../TMPResults/metrics/10x_v2_mousebrain/genescoring_RAGI_table.csv
