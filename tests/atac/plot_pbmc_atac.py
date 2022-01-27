import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from utils import tools as tl
import utils.tools as tt
from sklearn.preprocessing import LabelEncoder
from collections import Counter

base_path = '/Users/zhongyuanke/data/'
file_vcca = base_path + 'vcca/atac/pbmc_10k/vcca_temp.h5ad'

adata_vcca = sc.read_h5ad(file_vcca)

file_rna = '/Users/zhongyuanke/data/multimodal/atac_pbmc_10k/rna.h5ad'
rna = sc.read_h5ad(file_rna)
print(rna)
# sc.pp.filter_genes(rna, min_cells=2)
# sc.pp.log1p(rna)
# sc.pp.pca(rna)
# sc.pp.neighbors(rna)
# sc.tl.umap(rna)
# rna.obs['vcca_celltype_l2'] = adata_vcca.obs['vcca_celltype_l2']
# rna.write_h5ad(file_rna)
# sc.pl.umap(rna, color = 'vcca_celltype_l2')
# sc.pl.umap(rna, color='vcca_leiden')
# sc.pl.umap(rna, color = 'seurat_celltype')

# ---------------------------------------------------------------------------
file_atac = base_path + 'multimodal/atac_pbmc_10k/atac.h5ad'

atac = sc.read_h5ad(file_atac)
# print(atac)
# sc.pp.log1p(atac)
# sc.pp.filter_genes(atac, min_cells=2)
# sc.pp.pca(atac)
# sc.pp.neighbors(atac)
# sc.tl.umap(atac)
# atac.obs['vcca_celltype_l2'] = adata_vcca.obs['vcca_celltype_l2']
# atac.write_h5ad(file_atac)
print(atac)

# sc.pl.umap(atac, color='vcca_celltype_l2')
# sc.pl.umap(atac, color = 'vcca_leiden')

# ---------------------------------------------------------------------------
# wnn_umap_path='/Users/zhongyuanke/data/multimodal/atac_pbmc_10k/wnn_umap_filt.csv'
# wnn_umap = pd.read_csv(wnn_umap_path, index_col=0)
#
# wnn_umap = wnn_umap.values
# rna_emb = rna.obsm['X_umap']

# celltype = adata_vcca.obs['vcca_celltype_l2'].values
# encoder = LabelEncoder()
# label = encoder.fit_transform(celltype)


# c_map = 'tab20'
# size = 2
# xy_label_size = 10
# title_size = 14
# xy_label = 'UMAP_'

# fig = plt.figure(figsize=(40, 30))
# plt.subplots_adjust(wspace =0.3, hspace =0.3)
#
# ax = fig.add_subplot(111)
# plt.title('RNA', fontsize=title_size)
# ax.scatter(wnn_umap[:, 0], wnn_umap[:, 1], c=label, cmap=c_map, s=size,
#            linewidth=0)
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# plt.xlabel(xy_label + '1', size=xy_label_size)
# plt.ylabel(xy_label + '2', size=xy_label_size)
# plt.show()
# plt.close(fig)
# file_rna = '/Users/zhongyuanke/data/multimodal/atac_pbmc_10k/rna.h5ad'
# file_vcca = base_path + 'vcca/atac/pbmc_10k/vcca_temp.h5ad'
# adata_vcca = sc.read_h5ad(file_vcca)
#
# adata_vcca.obsm['X_umap']=wnn_umap
# sc.pl.umap(adata_vcca,color = ['vcca_celltype_l2'])

# ------------------------------------------


rna_emb = rna.obsm['X_umap']
atac_emb = atac.obsm['X_umap']
vcca_emb = adata_vcca.obsm['X_umap']

celltype = adata_vcca.obs['vcca_celltype_l2'].values
encoder = LabelEncoder()
label = encoder.fit_transform(celltype)
print(Counter(celltype))
print(Counter(label))
# label = list(map(int,adata_vcca.obs['leiden']))
# seurat_celltype_path = base_path + 'multimodal/atac_pbmc_10k/seurat_celltype_filt.csv'
#
# label = tl.get_label_by_txt(seurat_celltype_path)
# label = encoder.fit_transform(celltype)
# c_map = 'Dark2'
c_map = 'tab20b'

my_color_map = {
    0: 'green',
    1: 'royalblue',
    2: 'teal',
    3: 'limegreen',
    4: 'grey',
    5: 'gold',
    6: 'skyblue',
    7: 'darkslateblue',
    8: 'violet',
    9: 'blueviolet',
    10: 'sandybrown',
    11: 'olive',
    12: 'lightcoral',
    13: 'turquoise',
    14: 'slateblue',
    15: 'rosybrown',
    16: 'red',
    17: 'darkred',
    18: 'tan',
    19: 'saddlebrown',
    20: 'deepskyblue',
    21: 'cyan'
}

my_color_map = {
    0: 'green',
    1: 'royalblue',
    2: 'teal',
    3: 'limegreen',
    4: 'grey',
    5: 'skyblue',
    6: 'darkslateblue',
    7: 'violet',
    8: 'blueviolet',
    9: 'sandybrown',
    10: 'olive',
    11: 'lightcoral',
    12: 'turquoise',
    13: 'slateblue',
    14: 'rosybrown',
    15: 'red',
    16: 'darkred',
    17: 'tan',
    18: 'saddlebrown',
    19: 'deepskyblue',
    20: 'cyan'
}

label_color = []
for l in label:
    label_color.append(my_color_map[l])
size = 2
xy_label_size = 10
title_size = 14
xy_label = 'UMAP_'

fig = plt.figure(figsize=(50, 30))
plt.subplots_adjust(wspace =0.3, hspace =0.3)

ax = fig.add_subplot(131)
plt.title('RNA', fontsize=title_size)
ax.scatter(rna_emb[:, 0], rna_emb[:, 1], c=label_color, cmap=c_map, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel(xy_label + '1', size=xy_label_size)
plt.ylabel(xy_label + '2', size=xy_label_size)


ax = fig.add_subplot(132)
plt.title('ATAC', fontsize=title_size)
ax.scatter(atac_emb[:, 0], atac_emb[:, 1], c=label_color, cmap=c_map, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel(xy_label + '1', size=xy_label_size)
plt.ylabel(xy_label + '2', size=xy_label_size)

ax = fig.add_subplot(133)
plt.title('VCCA', fontsize=title_size)
ax.scatter(vcca_emb[:, 0], vcca_emb[:, 1], c=label_color, cmap=c_map, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel(xy_label + '1', size=xy_label_size)
plt.ylabel(xy_label + '2', size=xy_label_size)
plt.show()
plt.close(fig)