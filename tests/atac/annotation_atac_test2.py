import scanpy as sc
from utils import tools as tl
import umap
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.cluster import KMeans
from collections import Counter

base_path = '/Users/zhongyuanke/data/'


file_atac = '/Users/zhongyuanke/data/vcca/atac/pbmc_10k/pvcca_atac_01.h5ad'
file_activity = base_path + 'multimodal/atac_pbmc_10k/atac_filt.csv'
file_vcca = base_path + 'vcca/atac/pbmc_10k/pvcca_01.h5ad'
file_vcca2 = base_path + 'vcca/atac/pvcca_03.h5ad'
celltype_path = base_path + 'multimodal/atac_pbmc_10k/celltype_filt.csv'

file_rna = '/Users/zhongyuanke/data/multimodal/atac_pbmc_10k/rna.h5ad'
label = tl.get_label_by_txt(celltype_path)
# label = tl.text_label_to_number(label)

# file_rna = base_path + 'multimodal/protein/rna.csv'
# file_protein = base_path + 'multimodal/protein/protein.csv'
# file_vcca = base_path + 'vcca/pbmc_protein/pvcca_02.h5ad'
# file_vcca2 = base_path + 'vcca/pbmc_protein/pvcca_03.h5ad'

rna = sc.read_h5ad(file_rna)
adata = sc.read_h5ad(file_vcca)
# print(Counter(label))
# print(Counter(adata.obs['celltype']))
# adata.obs['seurat_celltype']=label
# adata.write_h5ad(file_vcca)

# sc.pp.log1p(adata)
# sc.tl.leiden(adata, resolution=1.4)
# sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
# sc.pl.rank_genes_groups(adata, n_genes=5, sharey=False,)
# --------------------------activity--------------------------------
adata_atac = sc.read_h5ad(file_atac)
# sc.pp.log1p(adata_atac)
# print(adata_atac)
# print(adata)
# adata_atac.obsm['X_vcca'] = adata.obsm['vcca']
# sc.pp.neighbors(adata_atac,use_rep='X_vcca')
# sc.tl.leiden(adata_atac,resolution=1.4)
#
# # adata_atac.write_h5ad('/Users/zhongyuanke/data/vcca/atac/pbmc_10k/pvcca_atac_01.h5ad')
# sc.tl.rank_genes_groups(adata_atac, 'leiden', method='t-test')
# sc.pl.rank_genes_groups(adata_atac, n_genes=5, sharey=False,)
# sc.tl.umap(adata_atac)
# sc.pl.umap(adata_atac, color=['leiden'], s=2, ncols=4)

# kmeans = KMeans(n_clusters=27, random_state=123, n_jobs=4).fit(adata.obsm['X_vcca'])
# kmeans_label = kmeans.labels_
# print(Counter(kmeans_label))
# kmeans_uns = pd.Series(kmeans_label, dtype='category')
# print(kmeans_uns)
# adata.obs['kmeans'] = kmeans_uns.values
# adata.uns['kmeans'] = kmeans_uns

print(adata)
# adata_test = sc.read_h5ad(file_vcca2)
# adata_protein = sc.read_csv(file_protein)

cluster2annotation_01_l1 = {
    '0': 'Classical Monocyte',
    '1': 'CD8 T-cell',
    '2': 'Classical Monocyte',
    '3': 'naive CD4 T-cell',
    '4': 'Memory CD4 T-cell',
    '5': 'Memory CD8 T-cell',
    '6': 'Classical Monocyte',
    '7': 'Treg',
    '8': 'naive CD4 T-cell',
    '9': 'Memory B cell',
    '10': 'gdT',
    '11': 'Non-classical Monocyte',
    '12': 'CD8 T-cell',
    '13': 'Naive B',
    '14': 'Classical Monocyte',
    '15': 'Treg',
    '16': 'CD8 T-cell',
    '17': 'MAIT',
    '18': 'Myeloid DC',
    '19': 'pDC',
    '20': 'Memory CD4 T-cell',
    '21': 'Memory B cell',
}
# cluster2annotation_kmeans_pvcca01 = {
#     '0': 'basophil',
#     '1': 'memory CD4 T-cell',
#     '2': 'basophil',
#     '3':'naive B-cell',
#     '4': 'classical monocyte',
#     '5': 'memory CD8 T-cell',
#     '6': 'neutrophil',
#     '7': 'memory B-cell',
#     '8': 'plasmacytoid DC',
#     '9': 'neutrophil',
#     '10': 'naive CD4 T-cell',
#     '11': 'memory CD4 T-cel',
#     '12': 'basophil',
#     '13': 'naive CD4 T-cell',
#     '14': 'memory CD8 T-cell',
#     '15': 'naive B-cell',
#     '16': 'basophil',
#     '17': 'naive CD4 T-cell',
#     '18': 'neutrophil',
#     '19': 'gdT',
#     '20': 'classical monocyte',
#     '21': 'neutrophil',
#     '22': 'basophil',
#     '23': 'NKT',
#     '24': 'memory B-cell',
#     '25': 'Monocyte',
#     '26': 'memory CD8 T-cell',
# }

cluster2annotation_01_l2 = {
    '0':'basophil',
    '1':'classical monocyte',
    '2':'naive B-cell',
    '3':'neutrophil',
    '4':'memory CD4 T-cell',
    '5':'naive CD4 T-cell',
    '6':'memory CD4 T-cell',
    '7':'memory CD8 T-cell',
    '8':'memory CD4 T-cell',
    '9':'Non-classical monocyte',
    '10':'NK-cell',
    '11':'classical monocyte',
    '12':'gdT-cell',
    '13':'memory B-cell 02',
    '14':'naive B-cell',
    '15':'Treg',
    '16':'classical monocyte-16',
    '17':'memory B-cell 01',
    '18':'gdT-cell',
    '19':'MAIT T-cell 01',
    '20':'neutrophil',
    '21':'MAIT T-cell-21',
    '22':'memory B-cell 01',
    '23':'naive CD4 T-cell',
    '24':'myeloid DC',
    '25':'naive CD4 T-cell',
    '26':'gdT-cell',
    '27':'plasmacytoid DC-27',
}

cluster2annotation_kmeans_pvcca01 = {
    0: 'basophil',
    1: 'memory CD4 T-cell',
    2: 'basophil',
    3: 'naive B-cell',
    4: 'classical monocyte',
    5: 'memory CD8 T-cell',
    6: 'neutrophil',
    7: 'memory B-cell',
    8: 'plasmacytoid DC',
    9: 'neutrophil',
    10: 'naive CD4 T-cell',
    11: 'memory CD4 T-cell',
    12: 'Myeloid DC',
    13: 'naive CD4 T-cell',
    14: 'memory CD8 T-cell',
    15: 'naive B-cell',
    16: 'Myeloid DC',
    17: 'naive CD4 T-cell',
    18: 'neutrophil',
    19: 'gdT',
    20: 'classical monocyte',
    21: 'neutrophil',
    22: 'basophil',
    23: 'Treg',
    24: 'memory B-cell',
    25: 'Monocyte',
    26: 'memory CD8 T-cell',
}

# sc.tl.leiden(adata)
# adata.write_h5ad(file_vcca)
# sc.pl.umap(adata, color=['seurat_celltype'], legend_loc='on data', frameon=False, s=5)

# sc.pl.umap(adata, color=['S100A8', 'CST3', 'CD8A', 'CD79A', 'CD74', 'MS4A1', 'LTB', 'FCER1G',
#                         'CD3E', 'CCL5', 'TYROBP', 'TCF7', 'LYZ', 'IL7R', 'HLA-DRA', 'LDHB'], s=2, ncols=4)
# sc.pp.neighbors(rna)
# sc.tl.umap(rna)
# sc.pl.umap(rna, color=['TRPC3'], s=6, ncols=4)
# sc.pl.umap(adata, color=['TMEM74'], s=8, ncols=4)
# sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
# sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False)

# adata.obs['wnn_label'] = label
# adata.obs['vcca_celltype'] = adata.obs['kmeans'].map(cluster2annotation_kmeans_pvcca01).astype('category')


adata.obs['vcca_celltype_l2'] = adata.obs['leiden'].map(cluster2annotation_01_l2).astype('category')
sc.pl.umap(adata, color='vcca_celltype_l2', legend_loc='on data', frameon=False, s=8)
# adata.obs['vcca_celltype_l2'] = adata.obs['leiden'].map(cluster2annotation_01_l2).astype('category')
# rna.obs['vcca_celltype_l2'] = adata.obs['vcca_celltype']
# rna.write_h5ad(file_rna)
# adata.write_h5ad(file_vcca)
adata_atac.obs['vcca_celltype_l2'] = adata_atac.obs['leiden'].map(cluster2annotation_01_l2).astype('category')
# sc.pl.umap(adata, color='vcca_celltype', legend_loc='on data', frameon=False, s=8)
adata_atac.write_h5ad(file_atac)
# print(adata)
# adata_test.obs['celltype'] = adata.obs['celltype']
# sc.pl.umap(adata_test, color='celltype', legend_loc='on data', frameon=False, s=10)
#
marker_map_rna = {
    'basophil':['RAB31'],
    'pDC': ['PLD4','TXNDC5'],
    'neutrophil': ['ARHGAP26','PLXDC2','CSF3R','SULF2','MCTP2'],
    'Non-classical Monocyte':['FCER1G','LYPD2'],
    'classical Monocyte':  ['FCN1', 'S1PR3', 'RXFP2','VCAN','DMXL2'],
    'gdT': ['CST7','NKG7',],
    'Treg':  ['IL7R', 'IL32', 'STAM'],
    'naive CD8 T-cell': 'CD8B',
    'naive CD4 T-cell': ['BCL11B', 'LEF1'],
    'Memory CD4 T-cell': 'INPP4B',
    'naive B': ['BACH2','SYN3','CD74','CD8A'],
    'NK T-cell': ['GNLY'],
    'MAIT': ['ME1','SLC4A10'],
    'memory B-cell': ['PPP1R37','BANK1','CD79A'],
    'memory CD8 T-cell':['CCL5','CD8A',],
    'Myeloid DC':['TPM2'],
}
marker_map_atac = {
    'basophil': ['RAB31'],
    'pDC': ['PLD4','TXNDC5','SOX4','ZFAT'],
    'neutrophil': ['MCTP2','PLXDC2','CSF3R','SULF2','MCTP2'],
    'Non-classical Monocyte':['FCER1G', 'AOAH', 'CABP4', 'LYPD2'],
    'classical Monocyte':  ['FCN1', 'S1PR3', 'RXFP2','VCAN','DMXL2'],
    'gdT': ['CST7','NKG7','JAKMIP1'],
    'Treg':  ['IL7R', 'IL32', 'STAM'],
    'CD8 T-cell': 'CD8B',
    'naive CD4 T-cell': ['BCL11B', 'LEF1'],
    'Memory CD4 T-cell': 'INPP4B',
    'naive B': ['BACH2','SYN3','CD74'],
    'NK T-cell': ['CACNA1C','GNLY'],
    'MAIT': ['ME1','SLC4A10','ZBTB16'],
    'memory B-cell': ['PPP1R37','BANK1','CD79A','PLEKHG1'],
    'memory CD8 T-cell':['CCL5','CD8A'],
    'Myeloid DC':['TPM2'],
}
# sc.pl.matrixplot(adata, marker_map, 'leiden', dendrogram=True, cmap='Blues',)

adata.layers['scaled'] = sc.pp.scale(adata, copy=True).X
sc.pl.matrixplot(adata, marker_map_rna, 'vcca_celltype_l2', dendrogram=True,
                 layer='scaled', vmin=-2, vmax=2, cmap='RdBu_r')

adata_atac.layers['scaled'] = sc.pp.scale(adata_atac, copy=True).X
sc.pl.matrixplot(adata_atac, marker_map_atac, 'vcca_celltype_l2', dendrogram=True,
                 layer='scaled', vmin=-2, vmax=2, cmap='RdBu_r')

# adata_atac.layers['scaled'] = sc.pp.scale(adata_atac, copy=True).X
# sc.pl.matrixplot(adata_atac, marker_map, 'vcca_celltype', dendrogram=False,
#                  layer='scaled', vmin=-2, vmax=2, cmap='RdBu_r')
# sc.pl.umap(adata, color=['SORCS1'],
#            s=2, ncols=4) # iMonocyte

# sc.pl.umap(adata, color=['CPA3'],
#            s=2, ncols=4) # balas
# sc.pl.umap(adata, color=['TPM2'],
#            s=2, ncols=4) # myeloid
# sc.pl.umap(adata, color=['THBS4'],
#            s=2, ncols=4) # Eosinophils

# sc.pl.dotplot(adata, marker_map,
#               groupby='wnn_label', dendrogram=True)

# ax = sc.pl.heatmap(adata, marker_map, groupby='celltype', cmap='viridis', dendrogram=True)
# sc.tl.rank_genes_groups(adata, 'celltype', method='t-test')
# sc.pl.rank_genes_groups(adata, n_genes=5, sharey=False)
# sc.pl.umap(adata, color=['MT-ATP6', 'CD74', 'LYZ', 'IL32', 'CD3D', 'CST7', 'LEF1', 'FCER1G', 'PLD4', 'FTL'],
#            s=2, ncols=4)

# sc.pl.umap(adata, color=['FCN1', 'CSF3R','ARHGAP25','ARHGAP26','SAMHD1','TYMP'],
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

