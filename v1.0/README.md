# GAGAM v1.0

GAGAM is Genomic Annotated Gene Activity Matrix, a model-driven approach to interpret scATAC-seq accessibility data, to obtain information about gene activity. It is a fundamental step to correlate scATAC-seq with ScRNA-seq. In other words, to link accessibility and expression of genes.

GAGAM v1.0 was used in the experimental setup published in:
Martini, L., Bardini, R., Savino, A., Di Carlo, S. (2022). GAGAM: A Genomic Annotation-Based Enrichment of scATAC-seq Data for Gene Activity Matrix. In: Rojas, I., Valenzuela, O., Rojas, F., Herrera, L.J., Ortuño, F. (eds) Bioinformatics and Biomedical Engineering. IWBBIO 2022. Lecture Notes in Computer Science(), vol 13347. Springer, Cham. https://doi.org/10.1007/978-3-031-07802-6_2

This README file lists all available scripts and the full workflow to follow to reproduce the experiment reported in the paper. 

Shey are divided into:

1. R scripts
1. Python scripts

## R scripts

1. **10x_v1_PBMChuman_processing.R**: contains the sripts to calculate GAGAM, Cicero and GeneScoring, for all the human PBMC datasets.

2. **buenrostro_processing.R** : contains the sripts to calculate GAGAM, Cicero and GeneScoring, for all the Buenrostro2018 dataset, containing bone marrow cells.

3. **snare_new_gam_processing.R**: contains the scripts to calculate GAGAM, Cicero and GeneScoring, for all the datasets containing mouse brain cells.

## Python scripts

Python scripts use the python3 language, and include:

1. **classifyPeaks.py**: labels peaks from a peak file using cCRECombined genomic annotation tracks;
2. **computeRAGI.py**: computes the RAGI scores for a GAM versus another classification, leveraging a list of housekeeping genes and a list of marker genes.

# Procedure

## Phase 0: Data preparation

For each dataset, download the data from the sources

- 10X V1.0.1 PBMC: https://www.10xgenomics.com/resources/datasets/5-k-peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-1-standard-1-0-1
- 10X V2.0.0 PBMC: https://www.10xgenomics.com/resources/datasets/5-k-peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-next-gem-v-1-1-1-1-standard-2-0-0
- Buenrostro2018: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96769
- 10X V1.1.0 Brain: https://www.10xgenomics.com/resources/datasets/fresh-cortex-from-adult-mouse-brain-p-50-1-standard-1-1-0
- SNARE: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126074

Organize the data following in propoided folder structure in the TMPDATA folder

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
- 10X_V1_Brain
  - atac_v1_adult_brain_fresh_5k_peaks.bed
  - atac_v1_adult_brain_fresh_5k_filtered_peak_bc_matrix.h5
- SNARE
  - GSE126074_AdBrainCortex_SNAREseq_chromatin.barcodes.tsv
  - GSE126074_AdBrainCortex_SNAREseq_chromatin.counts.mtx
  - GSE126074_AdBrainCortex_SNAREseq_chromatin.peaks.tsv

## Phase 1: Peak labeling

The first step is to label the peaks with Genomic annotation tracks.

## Input

- File containig peaks for each dataset
- Species appropriate cCRE Combined annotation tracks

## Script

- **classifyPeaks.py**

## Syntax

Execute the script from the v1.0 folder with the syntax:
```
##10X V1.0.1 PBMC peak labeling
python3 classifyPeaks.py '../TMPDATA/10x_V1_PBMChuman/peaks.bed' '../DATA/Genomic_annotations/encodeCcreCombined_human.bb'

##10X V2.0.0 PBMC peak labeling
python3 classifyPeaks.py '../TMPDATA/10X_V2_PBMC/peaks.bed' '../DATA/Genomic_annotations/encodeCcreCombined_human.bb'

##Buenrostro2018 peak labeling
python3 classifyPeaks.py '../TMPDATA/buenrostro/GSE96769_PeakFile_20160207.csv.bed' '../DATA/Genomic_annotations/encodeCcreCombined_human.bb'

##10X V1.1.0 Brain peak labeling
python3 classifyPeaks.py '../TMPDATA/10X_V1_Brain/atac_v1_adult_brain_fresh_5k_peaks.bed' '../DATA/Genomic_annotations/encodeCcreCombined_mouse.bb'

##SNARE peak labeling
python3 classifyPeaks.py '../TMPDATA/SNARE/GSE126074_AdBrainCortex_SNAREseq_chromatin.peaks.tsv' '../DATA/Genomic_annotations/encodeCcreCombined_mouse.bb'

```

## Output

The script will output a file containing the labeled peaks for each dataset in the respective folder:

- ../TMPResults/labeled_peaks/10x_V1_PBMChuman/cCRE_labeled_peaks.tsv
- ../TMPResults/labeled_peaks/10X_V2_PBMC/cCRE_labeled_peaks.tsv
- ../TMPResults/labeled_peaks/buenrostro/cCRE_labeled_peaks.tsv
- ../TMPResults/labeled_peaks/10X_V1_Brain/cCRE_labeled_peaks.tsv
- ../TMPResults/labeled_peaks/SNARE/cCRE_labeled_peaks.tsv

## Phase 2: GAMs construction and clustering analysis

The second phase consists in scATAC-seq data standard processing, construction of each of the three Gene Activity Matrices, and clustering classification of all of them.

### Input

Path for 10X V1.0.1 PBMC inputs:
  - ../TMPResults/leabeled_peaks/10x_V1_PBMChuman/encodeCcreCombined_hg38_ucscLabel_classifiedPeaks.csv
  - ../TMPDATA/10x_V1_PBMChuman/
  - ../DATA/Genomes/hg38/GCF_000001405.39_GRCh38.p13_genomic.gtf

Path for 10X V2.0.0 PBMC inputs:
  - ../TMPResults/labeled_peaks/10X_V2_PBMC/encodeCcreCombined_hg38_ucscLabel_classifiedPeaks.csv
  - ../TMPDATA/10X_V2_PBMC/
  - ../DATA/Genomes/hg38/GCF_000001405.39_GRCh38.p13_genomic.gtf

Path for Buenrostro2018 inputs:
  - ../TMPResults/labeled_peaks/buenrostro/encodeCcreCombined_hg38_ucscLabel_classifiedPeaks.csv
  - ../TMPDATA/buenrostro/
  - ../DATA/Genomes/hg38/GCF_000001405.39_GRCh38.p13_genomic.gtf

Path for 10X V1.1.0 Brain inputs:
  - ../TMPResults/labeled_peaks/10X_V1_Brain/encodeCcreCombined_ucscLabel_classifiedPeaks.csv
  - ../TMPDATA/10X_V1_Brain/
  - ../DATA/Genomes/mm10/GCF_000001635.26_GRCm38.p6_genomic.gtf

Path for SNARE inputs:
  - ../TMPResults/labeled_peaks/SNARE/encodeCcreCombined_ucscLabel_classifiedPeaks.csv
  - ../TMPDATA/SNARE/
  - ../DATA/Genomes/mm10/GCF_000001635.26_GRCm38.p6_genomic.gtf

### Script

- **10x_v1_PBMChuman_processing.R** 
- **10X_V2_PBMC_processing.R**
- **buenrostro_processing.R**
- **10X_V1_Brain.R**
- **snare_new_gam_processing.R**


### Syntax

Execute the script from the v1.0 folder with the syntax:
```
##10X V1.0.1 PBMC peak labeling
Rscript 10x_v1_PBMChuman_processing.R

##10X V2.0.0 PBMC peak labeling
Rscript 10X_V2_PBMC_processing.R

##Buenrostro2018 peak labeling
Rscript buenrostro_processing.R 

##10X V1.1.0 Brain peak labeling
Rscript 10X_V1_Brain.R 

##SNARE peak labeling
Rscript snare_new_gam_processing.R
```

### Output

The script execution outputs the three GAMs, their clustering classifications, scATAC-seq clustering classification and ARI and AMI tables. 

Path for 10x_v1_PBMChuman_processing.R outputs:
  - ../TMPResults/GAM/10x_V1_PBMChuman/gagam.tsv
  - ../TMPResults/GAM/10x_V1_PBMChuman/cicero.tsv
  - ../TMPResults/GAM/10x_V1_PBMChuman/genescoring.tsv
  - ../TMPResults/classifications/10x_V1_PBMChuman/gagam_classifications.tsv
  - ../TMPResults/classifications/10x_V1_PBMChuman/cicero_classifications.tsv
  - ../TMPResults/classifications/10x_V1_PBMChuman/genescoring_classifications.tsv
  - ../TMPResults/classifications/10x_V1_PBMChuman/atac_classification.tsv
  - ../TMPResults/metrics/10x_V1_PBMChuman/ARI_AMI_table.tsv

Path for 10X_V2_PBMC_processing.R outputs:
  - ../TMPResults/GAM/10X_V2_PBMC/gagam.tsv
  - ../TMPResults/GAM/10X_V2_PBMC/cicero.tsv
  - ../TMPResults/GAM/10X_V2_PBMC/genescoring.tsv
  - ../TMPResults/classifications/10X_V2_PBMC/gagam_classifications.tsv
  - ../TMPResults/classifications/10X_V2_PBMC/cicero_classifications.tsv
  - ../TMPResults/classifications/10X_V2_PBMC/genescoring_classifications.tsv
  - ../TMPResults/classifications/10X_V2_PBMC/atac_classification.tsv
  - ../TMPResults/metrics/10X_V2_PBMC/ARI_AMI_table.tsv

Path for buenrostro_processing.R  outputs:
  - ../TMPResults/GAM/buenrostro/gagam.tsv
  - ../TMPResults/GAM/buenrostro/cicero.tsv
  - ../TMPResults/GAM/buenrostro/genescoring.tsv
  - ../TMPResults/classifications/buenrostro/gagam_classifications.tsv
  - ../TMPResults/classifications/buenrostro/cicero_classifications.tsv
  - ../TMPResults/classifications/buenrostro/genescoring_classifications.tsv
  - ../TMPResults/classifications/buenrostro/atac_classification.tsv
  - ../TMPResults/metrics/buenrostro/ARI_AMI_table.tsv

Path for 10X_V1_Brain.R  outputs:
  - ../TMPResults/GAM/10X_V1_Brain/gagam.tsv
  - ../TMPResults/GAM/10X_V1_Brain/cicero.tsv
  - ../TMPResults/GAM/10X_V1_Brain/genescoring.tsv
  - ../TMPResults/classifications/10X_V1_Brain/gagam_classifications.tsv
  - ../TMPResults/classifications/10X_V1_Brain/cicero_classifications.tsv
  - ../TMPResults/classifications/10X_V1_Brain/genescoring_classifications.tsv
  - ../TMPResults/classifications/10X_V1_Brain/atac_classification.tsv
  - ../TMPResults/metrics/10X_V1_Brain/ARI_AMI_table.tsv

Path for snare_new_gam_processing.R outputs:
  - ../TMPResults/GAM/SNARE/gagam.tsv
  - ../TMPResults/GAM/SNARE/cicero.tsv
  - ../TMPResults/GAM/SNARE/genescoring.tsv
  - ../TMPResults/classifications/SNARE/gagam_classifications.tsv
  - ../TMPResults/classifications/SNARE/cicero_classifications.tsv
  - ../TMPResults/classifications/SNARE/genescoring_classifications.tsv
  - ../TMPResults/classifications/SNARE/atac_classification.tsv
  - ../TMPResults/metrics/SNARE/ARI_AMI_table.tsv

## Phase 4: Metrics Calculation

The forth phase is the calculation of RAGI scores between each GAM and all the classifications for each dataset. 

### Input

For each dataset, the script takes as inputs the file for the GAM under analysis, the path of the folder containing all classification results to perform the comparison, and paths to species- and tissue-specific housekeeping and marker genes list files.

### Script

- **computeRAGI.py**

### Syntax

Compute RAGI for each GAM over classification results, for each dataset.

```
#### 10X V1.0.1 PBMC ######
#GAGAM v1.0
python3 computeRAGI.py '../TMPResults/GAM/10x_V1_PBMChuman/gagam.tsv' '../TMPResults/classifications/10x_V1_PBMChuman/' '../DATA/IWBBIO_2022/Housekeeping_genes/Housekeeping_genes_human.txt' '../DATA/IWBBIO_2022/Marker_genes/human_PBMC_gene_markers.csv'
#Cicero GAM
python3 computeRAGI.py '../TMPResults/GAM/10x_V1_PBMChuman/cicero.tsv' '../TMPResults/classifications/10x_V1_PBMChuman/' '../DATA/IWBBIO_2022/Housekeeping_genes/Housekeeping_genes_human.txt' '../DATA/IWBBIO_2022/Marker_genes/human_PBMC_gene_markers.csv'
#GeneScoring GAM
python3 computeRAGI.py '../TMPResults/GAM/10x_V1_PBMChuman/genescoring.tsv' '../TMPResults/classifications/10x_V1_PBMChuman/' '../DATA/IWBBIO_2022/Housekeeping_genes/Housekeeping_genes_human.txt' '../DATA/IWBBIO_2022/Marker_genes/human_PBMC_gene_markers.csv'
```


```
#### 10X V2.0.0 PBMC ######
#GAGAM v1.0
python3 computeRAGI.py '../TMPResults/GAM/10X_V2_PBMC/gagam.tsv' '../TMPResults/classifications/10X_V2_PBMC/' '../DATA/IWBBIO_2022/Housekeeping_genes/Housekeeping_genes_human.txt' '../DATA/IWBBIO_2022/Marker_genes/human_PBMC_gene_markers.csv'
#Cicero GAM
python3 computeRAGI.py '../TMPResults/GAM/10X_V2_PBMC/cicero.tsv' '../TMPResults/classifications/10X_V2_PBMC/' '../DATA/IWBBIO_2022/Housekeeping_genes/Housekeeping_genes_human.txt' '../DATA/IWBBIO_2022/Marker_genes/human_PBMC_gene_markers.csv'
#GeneScoring GAM
python3 computeRAGI.py '../TMPResults/GAM/10X_V2_PBMC/genescoring.tsv' '../TMPResults/classifications/10X_V2_PBMC/' '../DATA/IWBBIO_2022/Housekeeping_genes/Housekeeping_genes_human.txt' '../DATA/IWBBIO_2022/Marker_genes/human_PBMC_gene_markers.csv'
```
```
#### Buenrostro2018 ######
#GAGAM v1.0
python3 computeRAGI.py '../TMPResults/GAM/buenrostro/gagam.tsv' '../TMPResults/classifications/buenrostro/' '../DATA/IWBBIO_2022/Housekeeping_genes/Housekeeping_genes_human.txt' '../DATA/IWBBIO_2022/Marker_genes/human_PBMC_and_bone_marrow_gene_markers.csv'
#Cicero GAM
python3 computeRAGI.py '../TMPResults/GAM/buenrostro/cicero.tsv' '../TMPResults/classifications/buenrostro/' '../DATA/IWBBIO_2022/Housekeeping_genes/Housekeeping_genes_human.txt' '../DATA/IWBBIO_2022/Marker_genes/human_PBMC_and_bone_marrow_gene_markers.csv'
#GeneScoring GAM
python3 computeRAGI.py '../TMPResults/GAM/buenrostro/genescoring.tsv' '../TMPResults/classifications/buenrostro/' '../DATA/IWBBIO_2022/Housekeeping_genes/Housekeeping_genes_human.txt' '../DATA/IWBBIO_2022/Marker_genes/human_PBMC_and_bone_marrow_gene_markers.csv'
```
```
#### 10X V1.1.0 Brain ######
#GAGAM v1.0
python3 computeRAGI.py '../TMPResults/GAM/10X_V1_Brain/gagam.tsv' '../TMPResults/classifications/10X_V1_Brain/' '../DATA/IWBBIO_2022/Housekeeping_genes/Housekeeping_genes_mouse.txt' '../DATA/IWBBIO_2022/Marker_genes/mouse_brain_gene_markers.csv'
#Cicero GAM
python3 computeRAGI.py '../TMPResults/GAM/10X_V1_Brain/cicero.tsv' '../TMPResults/classifications/10X_V1_Brain/' '../DATA/IWBBIO_2022/Housekeeping_genes/Housekeeping_genes_mouse.txt' '../DATA/IWBBIO_2022/Marker_genes/mouse_brain_gene_markers.csv'
#GeneScoring GAM
python3 computeRAGI.py '../TMPResults/GAM/10X_V1_Brain/genescoring.tsv' '../TMPResults/classifications/10X_V1_Brain/' '../DATA/IWBBIO_2022/Housekeeping_genes/Housekeeping_genes_mouse.txt' '../DATA/IWBBIO_2022/Marker_genes/mouse_brain_gene_markers.csv'
```
```
#### SNARE ######
#GAGAM v1.0
python3 computeRAGI.py '../TMPResults/GAM/SNARE/gagam.tsv' '../TMPResults/classifications/SNARE/' '../DATA/IWBBIO_2022/Housekeeping_genes/Housekeeping_genes_mouse.txt' '../DATA/IWBBIO_2022/Marker_genes/mouse_brain_gene_markers.csv'
#Cicero GAM
python3 computeRAGI.py '../TMPResults/GAM/SNARE/cicero.tsv' '../TMPResults/classifications/SNARE/' '../DATA/IWBBIO_2022/Housekeeping_genes/Housekeeping_genes_mouse.txt' '../DATA/IWBBIO_2022/Marker_genes/mouse_brain_gene_markers.csv'
#GeneScoring GAM
python3 computeRAGI.py '../TMPResults/GAM/SNARE/genescoring.tsv' '../TMPResults/classifications/SNARE/' '../DATA/IWBBIO_2022/Housekeeping_genes/Housekeeping_genes_mouse.txt' '../DATA/IWBBIO_2022/Marker_genes/mouse_brain_gene_markers.csv'
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

Path for 10X_V1_Brain  outputs:
  - ../TMPResults/metrics/10X_V1_Brain/gagam_RAGI_table.csv
  - ../TMPResults/metrics/10X_V1_Brain/cicero_RAGI_table.csv
  - ../TMPResults/metrics/10X_V1_Brain/genescoring_RAGI_table.csv

Path for SNARE outputs:
  - ../TMPResults/metrics/SNARE/gagam_RAGI_table.csv
  - ../TMPResults/metrics/SNARE/cicero_RAGI_table.csv
  - ../TMPResults/metrics/SNARE/genescoring_RAGI_table.csv
