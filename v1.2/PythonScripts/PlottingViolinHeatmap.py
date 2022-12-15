import sys
import scanpy as sc
import pandas as pd


def AEconfrontationPlot(A_AnnData, E_AnnData, gene_path, gene_10_path):

    adata = sc.read_h5ad(A_AnnData)
    print(adata)
    adata_RNA = sc.read_h5ad(E_AnnData)
    print(adata_RNA)
    adata_RNA.obs['cluster'] =adata_RNA.obs["activity"].copy()

    csv = pd.read_csv(gene_path)
    csv = csv['x']
    genes = list(csv)

    sc.pl.stacked_violin(adata, genes, groupby='cluster', swap_axes=True, cmap='viridis', title="ACTIVITY Violin Plot of marker genes", save="../TMPResults/IMAGES/Activity_violin_plot.png")
    sc.pl.stacked_violin(adata_RNA, genes, groupby='cluster',swap_axes=True, cmap='viridis', title="EXPRESSION Violin Plot of marker genes", save="../TMPResults/IMAGES/Expression_violin_plot.png")


    csv_10 = pd.read_csv(gene_10_path)
    csv_10 = csv_10['x']
    genes_10 = list(csv_10)

    sc.pl.heatmap(adata, genes_10, groupby='cluster',show_gene_labels=True, cmap='viridis', save="../TMPResults/IMAGES/Activity_heatmap.png")
    sc.pl.heatmap(adata_RNA, genes_10, groupby='cluster',show_gene_labels=True, cmap='viridis', save="../TMPResults/IMAGES/Expression_heatmap.png")


if __name__ == '__main__':

    # sys.argv[1] -> Activity h5ad file path
    # sys.argv[2] -> Expression h5ad file path
    # sys.argv[3] -> Gene list for violin plot file path
    # sys.argv[4] -> Gene list for heatmap plot h5ad file path

    Activity_AD = sys.argv[1]
    Expression_AD = sys.argv[2]
    genes_1 = sys.argv[3]
    genes_10 = sys.argv[4]

    AEconfrontationPlot(Activity_AD, Expression_AD, genes_1, genes_10)



