# GAGAM

GAGAM is Genomic Annotated Gene Activity Matrix, a model-driven approach to interpret scATAC-seq accessibility data, to obtain information about gene activity. It is a fundamental step to correlate scATAC-seq with ScRNA-seq. In other words, to link accessibility and expression of genes.

## Release notes

GAGAM v1.0: 

- First release of GAGAM.


GAGAM v1.2:
- Promoter peaks detection precision increased
- Changed Intragenic Matrix contribution. Now it considers only exons from the gene bodies. Name changed to Exon Matrix.
- Normalization of the final matrix.
- Added comparative analysis with scRNA-seq data.

## How to cite

- GAGAM v1.0: Martini, L., Bardini, R., Savino, A., Di Carlo, S. (2022). GAGAM: A Genomic Annotation-Based Enrichment of scATAC-seq Data for Gene Activity Matrix. In: Rojas, I., Valenzuela, O., Rojas, F., Herrera, L.J., Ortuño, F. (eds) Bioinformatics and Biomedical Engineering. IWBBIO 2022. Lecture Notes in Computer Science(), vol 13347. Springer, Cham. https://doi.org/10.1007/978-3-031-07802-6_2.

- GAGAM v1.2: Martini, L., Bardini, R., Savino, A., Di Carlo, S. (2022). GGAGAM v1.2: an improvement on peak labeling and Genomic Annotated Gene Activity Matrix construction. Submitted for review for MDPI Genes journal.

## Repo structure and organization

~~~~
├── README.md
├── DATA                                                   //Data files supporting GAGAM construction and metrics calculations
│   ├── GENE_2022                                          //Data files for the MDPI Genes 2022 publication
│   │   ├── Housekeeping_genes                             //Housekeeping genes lists for human and mouse for RAGI computation
│   │   └── Marker_genes                                   //Marker genes lists for human and mouse for RAGI computation
│   ├── Genomes                                            //Genomic assemblies data for GAGAM construction
│   ├── Genomic_annotations                                //cCRECombined genomic annotation tracks for human and mouse species for peak labeling
│   └── IWBBIO_2022                                        //Data files for the IWBBIO 2022 publication
│       ├── Housekeeping_genes                             //Housekeeping genes lists for human and mouse for RAGI computation
│       └── Marker_genes                                   //Marker genes lists for human and mouse for RAGI computation
├── RESULTS                                                //Published results
│   ├── GENE_2022                                          //Results shown in the MDPI Genes 2022 publication
│   │   ├── Clustering classification results              //clustering classification results for all the GAM considered
│   │   ├── IMAGES                                         //Images of the results in the publication
│   │   ├── Labeled_peaks                                  //Labeled peaks files from peak labeling
│   │   ├── RAGI                                           //RAGI scores from RAGI computation
│   │   └── Datasets metrics and histograms.xlsx           //Support for metrics scores visualization                          
│   └── IWBBIO_2022                                        //Results shown in the IWBBIO 2022 publication
│       ├── Labeled_peaks                                  //Labeled peaks files from peak labeling
│       ├── RAGI                                           //RAGI scores from RAGI computation
│       └── iwbbio_metrics_graphs.xlsx                     //Support for metrics scores visualization                          
├── v1.0                                                   //Scripts supporting the construction and analysis of GAGAM v1.0
│   ├── PythonScripts                                      //Python scripts
│   │   ├── classifyPeaks.py                               //Python script to label peaks with cCRE genomic annotations
│   │   └── computeRAGI.py                                 //Python script to compute RAGI scores of a GAM compared to a set of classifications
│   └── Rscripts                                           //R scripts
│   │   ├── 10x_v1_PBMChuman_processing.R                  //R script to preprocess PBMC datasets analysed in the IWBBIO 2022 publication
│   │   ├── buenrostro_processing.R                        //R script to preprocess the bone marrow dataset analysed in the IWBBIO 2022 publication
│   │   └── snare_new_gam_processing.R                     //R script to preprocess the SNARE dataset analysed in the IWBBIO 2022 publication
│   └── README.md                                          
└── v1.2
    ├── PythonScripts                                      //Python scripts
    │   ├── classifyPeaks.py                               //Python script to label peaks with cCRE genomic annotations
    │   └── computeRAGI.py                                 //Python script to compute RAGI scores of a GAM compared to a set of classifications
    └── Rscripts                                           //R scripts
    │   ├── 10x_v1_PBMChuman_processing.R                  //R script to process 10X V1.0.1 PBMC datasets analysed in the MDPI Genes 2022 publication
    │   ├── 10X_V2_PBMC_processing.R                       //R script to process 10X V2.0.0 PBMC dataset analysed in the MDPI Genes 2022 publication
    │   ├── buenrostro_processing.R                        //R script to process Buenrostro2018 dataset analysed in the MDPI Genes 2022 publication
    │   ├── 10x_multiome_processing.R                      //R script to process 10X Multiome Controller PBMC dataset analysed in the MDPI Genes 2022 publication
    │   ├── 10x_multiome_X_processing.R                    //R script to process 10X Multiome Chromium X PBMC analysed in the MDPI Genes 2022 publication
    │   ├── 10x_v2_mousebrain_processing.R                 //R script to process 10X V2.1.0 Brain dataset analysed in the MDPI Genes 2022 publication
    │   ├── GAGAM_compuation_hg38.R                        //R script to compute and preprocess GAGAM in case of data alligned to human hg38 genome assembly
    │   └── GAGAM_computation_mm10.R                       //R script to compute and preprocess GAGAM in case of data alligned to mouse mm10 genome assembly
    └── README.md                     
~~~~
