# GAGAM v1.0


GAGAM is Genomic Annotated Gene Activity Matrix, a model-driven approach to interpret scATAC-seq accessibility data, to obtain information about gene activity. It is a fundamental step to correlate scATAC-seq with ScRNA-seq. In other words, to link accessibility and expression of genes.

v1.0
Here there are all the scripts employed by GAGAM v1.0. They are divided in:

1. R scripts
1. Python scripts

## R scripts
There are all the R scripts useful to compute GAGAM v1.0, and specifically, to obtain the results shown in conference paper "GAGAM: A Genomic Annotation-Based Enrichment of scATAC-seq Data for Gene Activity Matrix" from IWBBIO 2022.

1. **10x_v1_PBMChuman_processing.R**: contains the sripts to calculate GAGAM, Cicero and GeneScoring, for all the human PBMC datasets.

2. **buenrostro_processing.R** : contains the sripts to calculate GAGAM, Cicero and GeneScoring, for all the Buenrostro2018 dataset, containing bone marrow cells.

3. **snare_new_gam_processing.R**: contains the scripts to calculate GAGAM, Cicero and GeneScoring, for all the datasets containing mouse brain cells.


These scripts takes in input:

1. The scATAC-seq data, specifying th path to all the elements (namely, Matrix, peaks and cells files)
2. The labeled_peaks.csv, output of **classifyPeaks.py**
3. The appropriate genome reference, available in the folder DATA/Genomes

The outputs are:

1. the GAM.tsv files, for GAGAM, Cicero, and GeneScoring
2. the classfication.tsv files, for all the classification to compare

Both are necessary for the **computeRAGI.py** scripts.

Execute, for example, the script from the v1.0 folder with the syntax:

```
Rscript Rscripts/10x_v1_PBMChuman_processing.R
```



## Python scripts

Python scripts use the python3 language, and include:

1. **classifyPeaks.py**: labels peaks from a peak file using cCRECombined genomic annotation tracks;
2. **computeRAGI.py**: computes the RAGI scores for a GAM versus another classification, leveraging a list of housekeeping genes and a list of marker genes.

The following sections explain how to launch the two scripts. 

Also, an example run for each script shows how to label peaks from the 10xv1_PBMC_Human dataset, and to compute RAGI scores between the relative GAGAM and the other classifications. Files supporting the example are:

1. The peak file for peak labeling, provided at: DATA/IWBBIO_2022/Peaks/10xv1_PBMC_Human.bed
2. The GAGAM file for RAGI computation, provided at: human_pbmc_example/10xv1_PBMC_Human_ga_gam1.tsv
3. The path to the folder containing classifications for RAGI computation, provided at human_pbmc_example/classifications/10xv1_PBMC_Human_classifications

### classifyPeaks.py

This script takes as arguments:

1. the filename of peaks to be classified (\<peaksFilename\>)
2. the filename of the genomic annotation track (\<annotationTrackFilename\>)

By default, the script uses 'uscsLabel' annotations to label peaks.

Execute the script from the v1.0 folder with the syntax:

```
python3 classifyPeaks.py <peaksFilename> <annotationTrackFilename>
```


### computeRAGI.py

This script takes as arguments:

1. the filename of the GAM to consider (\<GAMFilename\>)
2. the path of the folder containing classification results to consider (\<classificationsFolderPath\>)
3. the filename of the housekeeping genes list to consider (\<housekeepingFilename\>)
4. the filename of the marker genes list to consider (\<markerFilename\>)

Execute the script from the v1.0 folder with the syntax:

```
python3 computeRAGI.py <GAMFilename> <classificationsFolderPath> <housekeepingFilename> <markerFilename>
```

# Procedure

## Phase 1: Peak labeling

THe first step ois to label the peaks with Genomic annotation tracks.

**Input**
 
- peaks_File
- annotationTrack_File

**Script**

- **classifyPeaks.py**

**Syntax**

Execute the script from the v1.0 folder with the syntax:
```
python3 classifyPeaks.py <peaksFilename> <annotationTrackFilename>
```

**Output**

- labeled_peaks

## Phase 2: GAMs construction

The second phase consists in ScATAC-seq data standard processing, and construction of each of the three Gene Activity Matrices.

### Input

- labeled_peaks
- ScATAC-seq data
- Genomic_Annotation

### Script

- **10x_v1_PBMChuman_processing.R** 
- **buenrostro_processing.R**
- **snare_new_gam_processing.R**

### Syntax

Execute ny of the scripts from the v1.0/Rscripts folder with the syntax:
```
Rscript <script_processing.R> 
```

### Output

- GAGAM
- Cicero_GAM
- GeneScoring_GAM

## Phase 3: Classification

The third phase follows with the cells classification based on each GAm and ScATAC-seq data

### Input

- GAGAM
- Cicero_GAM
- GeneScoring_GAM
- ScATAC-seq data

### Script

- **10x_v1_PBMChuman_processing.R** 
- **buenrostro_processing.R**
- **snare_new_gam_processing.R**

### Syntax

Execute any of the scripts from the v1.0/Rscripts folder with the syntax:
```
Rscript <script_processing.R> 
```
### Output

- GAGAM_clustering
- Cicero_GAM_clustering
- GeneScoring_GAM_clustering
- ATAC_clustering


## Phase 4: Metrics Calculation

The forth phase is the calculation of the ARI, AMI, and RAGI between GAGAM_clustering and all the other classifications. 

### Input

- GAGAM_clustering
- Cicero_GAM_clustering
- GeneScoring_GAM_clustering
- ATAC_clustering

### Script

- **computeRAGI.py**
- **10x_v1_PBMChuman_processing.R** 
- **buenrostro_processing.R**
- **snare_new_gam_processing.R**


### Syntax

Execute any of the scripts from the v1.0/Rscripts folder with the syntax:
```
Rscript <script_processing.R> 
```

and

```
python3 computeRAGI.py <GAMFilename> <classificationsFolderPath> <housekeepingFilename> <markerFilename>
```

### Output

- ARI, AMI, RAGI metric tables