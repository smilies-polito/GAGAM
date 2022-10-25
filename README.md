# GAGAM

GAGAM is Genomic Annotated Gene Activity Matrix, a model-driven approach to interpret scATAC-seq accessibility data, to obtain information about gene activity. It is a fundamental step to correlate scATAC-seq with ScRNA-seq. In other words, to link accessibility and expression of genes.

## Release notes

GAGAM v1.0: 

- First release of GAGAM.

## How to cite

- GAGAM v1.0: Martini, L.; Bardini, R.; Savino, A.; Di Carlo, S. "GAGAM: A Genomic Annotation-Based Enrichment of scATAC-seq Data for Gene Activity Matrix". In Proceedings of the Bioinformatics and Biomedical Engineering;


## Repo structure and organization

~~~~
├── README.md
├── DATA                                                   //Data files supporting GAGAM construction and metrics calculations
│   ├── Genomes                                            //Genomic assemblies data for GAGAM construction
│   ├── Genomic_annotations                                //cCRECombined genomic annotation tracks for human and mouse species for peak labeling
│   └── IWBBIO_2022                                        //Data files for the IWBBIO 2022 publication
│       ├── Housekeeping_genes                             //Housekeeping genes lists for human and mouse for RAGI computation
│       └── Marker_genes                                   //Marker genes lists for human and mouse for RAGI computation
├── RESULTS                                                //Published results
│   └── IWBBIO_2022                                        //Results shown in the IWBBIO 2022 publication
│       ├── Labeled_peaks                                  //Labeled peaks files from peak labeling
│       ├── RAGI                                           //RAGI scores from RAGI computation
│       └── iwbbio_metrics_graphs.xlsx                     //Support for metrics scores visualization                          
└── v1.0                                                   //Scripts supporting the construction and analysis of GAGAM v1.0
    ├── PythonScripts                                      //Python scripts
    │   ├── classifyPeaks.py                               //Python script to label peaks with cCRE genomic annotations
    │   └── computeRAGI.py                                 //Python script to compute RAGI scores of a GAM compared to a set of clas│sifications
    └── Rscripts                                           //R scripts
    │   ├── 10x_v1_PBMChuman_processing.R                  //R script to preprocess PBMC datasets analysed in the IWBBIO 2022 publication
    │   ├── buenrostro_processing.R                        //R script to preprocess the bone marrow dataset analysed in the IWBBIO 2022 publication
    │   └── snare_new_gam_processing.R                     //R script to preprocess the SNARE dataset analysed in the IWBBIO 2022 publication
    └── README.md                                          

    
~~~~
