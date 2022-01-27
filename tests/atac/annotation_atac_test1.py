import scanpy as sc
from utils import tools as tl
import umap
import numpy as np
import matplotlib.pyplot as plt


base_path = '/Users/zhongyuanke/data/'


file_rna = base_path + 'multimodal/atac_pbmc_10k/rna_filt.csv'
file_protein = base_path + 'multimodal/atac_pbmc_10k/atac.h5ad'
file_vcca = base_path + 'vcca/atac/pbmc_10k/pvcca_01.h5ad'
file_vcca2 = base_path + 'vcca/atac/pvcca_03.h5ad'
celltype_path = base_path + 'multimodal/atac_pbmc_10k/celltype_filt.csv'

label = tl.get_label_by_txt(celltype_path)
# label = tl.text_label_to_number(label)

# file_rna = base_path + 'multimodal/protein/rna.csv'
# file_protein = base_path + 'multimodal/protein/protein.csv'
# file_vcca = base_path + 'vcca/pbmc_protein/pvcca_02.h5ad'
# file_vcca2 = base_path + 'vcca/pbmc_protein/pvcca_03.h5ad'


adata = sc.read_h5ad(file_vcca)
rna = sc.read_csv(file_rna)
print(adata)
# adata_test = sc.read_h5ad(file_vcca2)
# adata_protein = sc.read_csv(file_protein)

cluster2annotation_02_l2 = {
    '0': 'classical Monocyte',
    '1': 'CD8 T-cell',
    '2': 'Neutrophil',
    '3': 'naive CD4 T-cell',
    '4': 'Memory CD4 T',
    '5': 'Treg',
    '6': 'classical Monocyte',
    '7': 'Treg',
    '8': 'naive CD4 T-cell',
    '9': 'B cell',
    '10': 'gdT',
    '11': 'Non-classical Monocyte',
    '12': 'CD8 T-cell',
    '13': 'Naive B',
    '14': 'classical Monocyte',
    '15': 'Treg',
    '16': 'CD8 T-cell',
    '17': 'MAIT',
    '18': 'Neutrophil',
    '19': 'pDC',
    '20': 'Memory CD4 T-cell',
    '21': 'B cell',
}

# sc.tl.leiden(adata)
# adata.write_h5ad(file_vcca)
# sc.pl.umap(adata, color='leiden', legend_loc='on data', frameon=False, s=5)
# sc.pl.umap(adata, color=['S100A8', 'CST3', 'CD8A', 'CD79A', 'CD74', 'MS4A1', 'LTB', 'FCER1G',
#                         'CD3E', 'CCL5', 'TYROBP', 'TCF7', 'LYZ', 'IL7R', 'HLA-DRA', 'LDHB'], s=2, ncols=4)
# sc.pp.neighbors(rna)
# sc.tl.umap(rna)
# sc.pl.umap(rna, color=['TRPC3'], s=6, ncols=4)
# sc.pl.umap(adata, color=['TMEM74'], s=8, ncols=4)
# sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
# sc.pl.rank_genes_groups(adata, n_genes=19, sharey=False)
adata.obs['wnn_label'] = label
adata.obs['celltype'] = adata.obs['leiden'].map(cluster2annotation_02_l2).astype('category')
sc.pl.umap(adata, color='celltype', legend_loc='on data', frameon=False, s=5)
# adata.write_h5ad(file_vcca)
# print(adata)
# adata_test.obs['celltype'] = adata.obs['celltype']
# sc.pl.umap(adata_test, color='celltype', legend_loc='on data', frameon=False, s=10)
#
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
# sc.pl.matrixplot(adata, marker_map, 'leiden', dendrogram=True, cmap='Blues',)
adata.layers['scaled'] = sc.pp.scale(adata, copy=True).X
sc.pl.matrixplot(adata, marker_map, 'celltype', dendrogram=True,
                 layer='scaled', vmin=-2, vmax=2, cmap='RdBu_r')
# sc.pl.dotplot(adata, marker_map,
#               groupby='wnn_label', dendrogram=True)

# ax = sc.pl.heatmap(adata, marker_map, groupby='celltype', cmap='viridis', dendrogram=True)
# sc.tl.rank_genes_groups(adata, 'celltype', method='t-test')
# sc.pl.rank_genes_groups(adata, n_genes=5, sharey=False)
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

