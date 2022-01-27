import scanpy as sc
from utils import tools as tl
import umap
import numpy as np
import matplotlib.pyplot as plt


base_path = '/Users/zhongyuanke/data/'

file_rna = base_path + 'multimodal/atac_pbmc_10k/rna_filt.csv'
file_atac = base_path + 'multimodal/atac_pbmc_10k/atac.h5ad'
celltype_path = base_path + 'multimodal/atac_pbmc_10k/celltype_filt.csv'

# file_rna = base_path + 'multimodal/atac/rna.csv'
# file_protein = base_path + 'multimodal/atac/protein.csv'
file_vcca = base_path + 'vcca/atac/pbmc_10k/pvcca_01.h5ad'
file_vcca2 = base_path + 'vcca/atac/pvcca_03.h5ad'


# file_rna = base_path + 'multimodal/protein/rna.csv'
# file_protein = base_path + 'multimodal/protein/protein.csv'
# file_vcca = base_path + 'vcca/pbmc_protein/pvcca_02.h5ad'
# file_vcca2 = base_path + 'vcca/pbmc_protein/pvcca_03.h5ad'


adata = sc.read_h5ad(file_vcca)
# adata_test = sc.read_h5ad(file_vcca2)
# adata_protein = sc.read_csv(file_protein)

cluster2annotation_01_l1 = {
    '0': 'cDC',
    '1': 'T cell',
    '2': 'T cell',
    '3': 'T cell',
    '4': 'cDC',
    '5': 'T cell',
    '6': 'MAIT',
    '7': 'cDC',
    '8': 'MAIT',
    '9': 'T cell',
    '10': 'NK',
    '11': 'NK',
    '12': 'T cell',
    '13': 'NK2',
    '14': 'cDC',
    '15': 'NK2',
}
cluster2annotation_01_l2 = {
    '0': 'B cell',
    '1': 'Naive CD4+ T cell',
    '2': 'CD4+Cytotoxic T cell',
    '3': 'NK T cell',
    '4': 'B cell',
    '5': 'Treg',
    '6': 'NK cell',
    '7': 'B cell',
    '8': 'NK cell',
    '9': 'NK T cell',
    '10': 'NK cell',
    '11': 'NK T cell',
    '12': 'CD4+ Memory T cell',
    '13': 'B cell',
    '14': 'cDC(CD1C+)',
    '15': 'B cell',
}

cluster2annotation_02_l1 = {
    '0': 'Treg',
    '1': 'Monocyte',
    '2': 'B cell',
    '3': 'naive CD4 T-cell',
    '4': 'Monocyte',
    '5': 'gdT',
    '6': 'Monocyte',
    '7': 'Treg',
    '8': 'Monocyte',
    '9': 'gdT',
    '10': 'pDC',
    '11': 'B cell',
    '12': 'Treg',
    '13': 'CD8 T-cell',
    '14': 'myeloid',
    '15': 'gdT',
    '16': 'Treg',
    '17': 'NK T-cell'
}

cluster2annotation_02_l2 = {
    '0': 'Treg',
    '1': 'classical Monocyte',
    '2': 'B cell',
    '3': 'naive CD4 T-cell',
    '4': 'classical Monocyte',
    '5': 'gdT',
    '6': 'classical Monocyte',
    '7': 'Treg',
    '8': 'Non-classical Monocyte',
    '9': 'gdT',
    '10': 'pDC',
    '11': 'B cell',
    '12': 'Treg',
    '13': 'CD8 T-cell',
    '14': 'Non-classical Monocyte',
    '15': 'gdT',
    '16': 'Treg',
    '17': 'B cell'
}


sc.pl.umap(adata, color='louvain', legend_loc='on data', frameon=False, s=5)
sc.pl.umap(adata, color=['S100A8', 'CST3', 'CD8A', 'CD79A', 'CD74', 'MS4A1', 'LTB', 'FCER1G',
                        'CD3E', 'CCL5', 'TYROBP', 'TCF7', 'LYZ', 'IL7R', 'HLA-DRA', 'LDHB'], s=2, ncols=4)

# sc.pl.umap(adata, color=['LTB'], s=6, ncols=4)
# sc.tl.rank_genes_groups(adata, 'louvain', method='t-test')
# sc.pl.rank_genes_groups(adata, n_genes=19, sharey=False)

adata.obs['celltype'] = adata.obs['louvain'].map(cluster2annotation_02_l2).astype('category')
# sc.pl.umap(adata, color='celltype', legend_loc='on data', frameon=False, s=5)
# adata.write_h5ad(file_vcca)
# print(adata)
# adata_test.obs['celltype'] = adata.obs['celltype']
# sc.pl.umap(adata_test, color='celltype', legend_loc='on data', frameon=False, s=10)
# #
# marker_map = {
#     'pDC': 'PLD4',
#     'Non-classical Monocyte': 'FCER1G',
#     'classical Monocyte': 'FCN1',
#
#     'B cell': ['CD79A'],
#     'gdT': 'CST7',
#     'Treg':  ['IL7R', 'IL32', ],
#     'CD8 T-cell': 'CD8B',
#     'naive CD4 T-cell': 'FHIT',
#     # 'NK T-cell': ['CD8A', 'CD8B'],
# }

marker_map = {
    'pDC': ['PLD4'],
    'Non-classical Monocyte': ['FCER1G', 'AOAH', 'CABP4', 'LYPD2'],
    'classical Monocyte': ['FCN1', 'S1PR3', 'RXFP2'],
    'B cell': ['CD79A'],
    'gdT': 'CST7',
    'Treg':  ['IL7R', 'IL32', ],
    'CD8 T-cell': 'CD8B',
    'naive CD4 T-cell': 'FHIT',
    'Memory CD4 T-cell': 'INPP4B',
    'naive B': ['BACH2','SYN3'],
    'NK T-cell': ['CD8A', 'CD8B'],
    'MAIT': ['ME1','SLC4A10'],
    'memory B': 'PPP1R37',
}

# sc.pl.matrixplot(adata, marker_map, 'celltype', dendrogram=True, cmap='Blues',
#                  colorbar_title='column scaled\nexpression')
adata.layers['scaled'] = sc.pp.scale(adata, copy=True).X
sc.pl.matrixplot(adata, marker_map, 'celltype', dendrogram=True,
                 layer='scaled', vmin=-2, vmax=2, cmap='RdBu_r')
# sc.pl.dotplot(adata, marker_map,
#               groupby='celltype', dendrogram=True)

# ax = sc.pl.heatmap(adata, marker_map, groupby='celltype', cmap='viridis', dendrogram=True)
# sc.tl.rank_genes_groups(adata, 'celltype', method='t-test')
# sc.pl.rank_genes_groups(adata, n_genes=18, sharey=False)
# sc.pl.umap(adata, color=['CD79A', 'CD8B', 'LYZ', 'IL32', 'CD3D', 'CST7', 'LEF1', 'FCER1G', 'PLD4', 'FTL'],
#            s=2, ncols=4)


# adata_protein.obs['celltype'] = adata.obs['celltype']
# sc.pp.neighbors(adata_protein)
# sc.tl.umap(adata_protein)
# sc.pl.umap(adata_protein, color='celltype', legend_loc='on data', frameon=False, s=5)
#
# sc.pp.neighbors(adata)
# sc.tl.umap(adata)
# sc.pl.umap(adata, color='celltype', legend_loc='on data', frameon=False, s=5)
# sc.pl.umap(adata, color=['S100A8', 'CST3', 'CD8A', 'CD79A', 'CD74', 'MS4A1', 'LTB', 'FCER1G',
#                        'CD3E', 'CCL5', 'TYROBP', 'TCF7', 'LYZ', 'IL7R', 'HLA-DRA', 'LDHB'], s=2, ncols=4)

