import sys

import pandas as pd
import numpy as np
from scipy import stats
import os
import time

# The createBiomarkersList function returns a list of unique marker genes
# taking as input CellMarkerFile .csv file downloaded
# from the CellMarker 2.0 database(http://bio-bigdata.hrbmu.edu.cn/CellMarker/)
# performing a search based on "Species" and "Tissue type" relevant to the dataset
# under analysis, and considering ALL "Cell Names"
def createBiomarkersList(CellMarkerFile):

    # extracting data from CellMarker file
    biomarkers_df = pd.read_csv(CellMarkerFile)
    # filtering out markers from cancer samples
    biomarkers_df = biomarkers_df[biomarkers_df['Cancer']=='Normal']
    # considering only the column with lists of markers
    biomarkers_df = biomarkers_df['Cell Marker']

    # unpacking the obtained list of lists in a single, possibly redundant list
    biomarkers_str_list = np.array(list(biomarkers_df))
    biomarkers_list = np.array([])

    for i in biomarkers_str_list:
        biomarkers_list = np.append(i.split(','), biomarkers_list)

    # considering unique markers
    biomarkers_unique = np.unique(biomarkers_list)

    # return unique markers list
    return biomarkers_unique

# The createHousekeepingGenesList function returns a list of housekeeping genes
# taking as input a .tsv file containing a plain list of genes in column
# TODO: insert source of housekeeping genes used
def createHousekeepingGenesList(HousekeepingGenesFile):

    # extracting the list of housekeeping genes
    housekeepingGenes_df = pd.read_csv(HousekeepingGenesFile, delimiter='\t')
    # consider the first column, containing the genes list and make it a list
    #housekeepingGenes_df = housekeepingGenes_df[0]
    housekeepingGenes_list = list(housekeepingGenes_df)

    # return the list of housekeeping genes
    return housekeepingGenes_list


# The residual_average_gini_index computes RAGI taking as arguments:
# genes_scores -> the GAM
# folder_clusters -> the path to the folder containing classifications
# housekeeping_genes -> the housekeeping genes list
# marker_genes -> the marker genes list
# min_cells_per_cluster -> a parameter to filter out clusters smaller than n cells
def residual_average_gini_index(gene_scores,folder_clusters,
                                housekeeping_genes,marker_genes,
                                min_cells_per_cluster=10):

    # Creating a subset from the GAM including only the housekeeping genes and marker genes
    df_matrix_housekeeping = gene_scores.loc[gene_scores.index.intersection(housekeeping_genes),]
    df_matrix_marker = gene_scores.loc[gene_scores.index.intersection(marker_genes),]

    # The gini function computes the Gini score
    # given a list of values
    def gini(list_of_values):
        sorted_list = sorted(list_of_values)
        height, area = 0, 0
        for value in sorted_list:
            height += value
            area += height - value / 2.
            fair_area = height * len(list_of_values) / 2.
        return np.divide(np.subtract(fair_area, area), fair_area, out = None, where = fair_area != 0)

    # The calculate_gini function calculates Gini values for a gene
    def calculate_gini(df_matrix, gene_name, clustering_info):
        return gini(get_avg_per_cluster(df_matrix,gene_name,clustering_info,use_log2=False))

    # The calculate_gini_values function calculates Gini value for all the considered genes
    def calculate_gini_values(df_matrix, clustering_info):
        gini_values=[]
        for gene_name in df_matrix.index:
            gini_values.append(calculate_gini(df_matrix, gene_name,clustering_info))
        return gini_values

    # The score_clustering_solution function computes delta difference
    # of the average accessibility in marker vs housekeeping genes
    # and results of the Kolmogorov Smirnov test
    def score_clustering_solution(df_matrix_marker,df_matrix_housekeeping,clustering_info):
        gini_values_housekeeping=calculate_gini_values(df_matrix_housekeeping,clustering_info)
        gini_values_marker=calculate_gini_values(df_matrix_marker,clustering_info)
        statistic,p_value=stats.ks_2samp(gini_values_marker,gini_values_housekeeping)

        return  np.mean(gini_values_marker), np.mean(gini_values_housekeeping), np.mean(gini_values_marker)-np.mean(gini_values_housekeeping), statistic, p_value

    # The get_avg_per_cluster function computes the average accessibility value per cluster
    def get_avg_per_cluster(df_matrix, gene_name, clustering_info,use_log2=False):
        N_clusters=len(clustering_info.index.unique())
        avg_per_cluster=np.zeros(N_clusters)
        for idx,idx_cluster in enumerate(sorted(np.unique(clustering_info.index.unique()))):
            if use_log2:
                values_cluster=df_matrix.loc[gene_name,clustering_info.loc[idx_cluster,:].values.flatten()].apply(lambda x:np.log2(x+1))
            else:
                values_cluster=df_matrix.loc[gene_name,clustering_info.loc[idx_cluster,:].values.flatten()]

            avg_per_cluster[idx]=values_cluster.mean()
            if avg_per_cluster[idx]>0:
                  avg_per_cluster[idx]=avg_per_cluster[idx]

        return avg_per_cluster


    # Computing RAGI for the GAM versus all the other classifications

    # Initialize metrics dataframe
    df_metrics = pd.DataFrame(columns=['Method', 'Clustering', 'Gini_Marker_Genes', 'Gini_Housekeeping_Genes', 'Difference', 'KS_statistics', 'p-value'])

    # Iterate RAGI computation over the classifications at the folder_clusters path
    for clusters_filename in os.listdir(folder_clusters):

        method = '_'.join(clusters_filename.split('_')[:-1])
        print('Computing RAGI over ', method)

        # Extracting classification labels from the current classification
        df_clusters = pd.read_csv(os.path.join(folder_clusters, clusters_filename), sep = '\t', index_col = 0)

        for clustering_method in df_clusters.columns:

            clustering_info = pd.DataFrame(df_clusters[clustering_method])
            clustering_info['Barcode'] = clustering_info.index
            clustering_info = clustering_info.set_index(clustering_method)

            # Removing clusters with less than min_cells_per_cluster cells in the current classification
            cluster_sizes = pd.value_counts(clustering_info.index)
            clustering_info = clustering_info.loc[cluster_sizes[cluster_sizes > min_cells_per_cluster].index.values, :]

            # Compute RAGI for the current classification
            mean_gini_marker, mean_gini_housekeeping, mean_gini_difference, statistics, p_value = score_clustering_solution(df_matrix_marker, df_matrix_housekeeping, clustering_info)

            df_metrics = df_metrics.append({'Method': method, 'Clustering':clustering_method,
                               'Gini_Marker_Genes':mean_gini_marker, 'Gini_Housekeeping_Genes':mean_gini_housekeeping,
                               'Difference':mean_gini_difference, 'KS_statistics':statistics, 'p-value':p_value},
                              ignore_index=True)

    return df_metrics


if __name__ == '__main__':

    # sys.argv[1] -> <GAMFilename>
    # sys.argv[2] -> <classificationsFolderPath>
    # sys.argv[3] -> <housekeepingFilename>
    # sys.argv[4] -> <markerFilename>

    # path to GAM file
    gene_activity_matrix_path = '../human_pbmc_example/' + sys.argv[1] #project_name + '/computing_RAGI/gene_scores/'
    print ('Considering GAM file at: ' + gene_activity_matrix_path)


    # path to classifications folder
    clustering_outputs_folder = '../human_pbmc_example/classifications/' + sys.argv[2] + '/'  #project_name + '/computing_RAGI/clustering_outputs/'
    print ('Considering classifications at: ' + clustering_outputs_folder)

    # path to housekeeping genes file
    housekeeping_genes_path = '../DATA/IWBBIO_2022/Housekeeping_genes/' + sys.argv[3] #project_name + '/computing_RAGI/housekeeping_genes/'
    print ('Considering housekeeping genes file at: ' + housekeeping_genes_path)

    # path to housekeeping genes file
    marker_genes_path = '../DATA/IWBBIO_2022/Marker_genes/' + sys.argv[4] #project_name + '/computing_RAGI/marker_genes/'
    print ('Considering marker genes file at: ' + marker_genes_path)

    # path to RAGI computations results
    results_RAGI_path = '../human_pbmc_example/Performance_metrics/RAGI/' #project_name + '/computing_RAGI/RAGI/'
    print ('RAGI computation results will be saved at: ' + results_RAGI_path)

    # extracting housekeeping genes list
    # data from TODO: insert link
    housekeeping_genes_list = createHousekeepingGenesList(housekeeping_genes_path)# + '/housekeeping_genes.txt')
    print('Creating housekeeping genes list from: ' + housekeeping_genes_path)# + '/housekeeping_genes.txt')

    # extracting marker genes list
    # data from http://bio-bigdata.hrbmu.edu.cn/CellMarker
    marker_genes_list = createBiomarkersList(marker_genes_path)#'genes_journal22/computing_RAGI/references/human_PBMC_pure_pliner/human_PBMC_marker_genes_pure_pliner.csv')
    print('Creating marker genes list from: ' + marker_genes_path)#print('Considering marker genes from the list: ' + marker_genes_path + '/marker_genes.csv')

    # reading GAM file
    gene_scores = pd.read_csv(gene_activity_matrix_path, sep ='\t', index_col=0)
    print('Extracting gene activity scores from: ' + gene_activity_matrix_path)

    # computing RAGI for GAM file versus all clustering classification in clustering_output_folder
    df_metrics = residual_average_gini_index(gene_scores, clustering_outputs_folder, housekeeping_genes_list, marker_genes_list, min_cells_per_cluster=10)

    # save RAGI results file
    if os.path.isdir(results_RAGI_path):
        df_metrics.to_csv(results_RAGI_path + '/RAGI_results.csv') #+ os.path.splitext(gene_activity_matrix_path)[0] + '.csv')
    else:
        os.makedirs(results_RAGI_path)
        df_metrics.to_csv(results_RAGI_path + '/RAGI_results.csv')#+ os.path.splitext(gene_activity_matrix_path)[0] + '.csv')

    print('Finished RAGI computations for ' + gene_activity_matrix_path)
