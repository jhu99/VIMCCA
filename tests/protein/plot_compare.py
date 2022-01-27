import scanpy as sc
import matplotlib.pyplot as plt
from utils import tools as tl
import numpy as np
import argparse
from sklearn.preprocessing import normalize
import umap
from sklearn.manifold import TSNE
from sklearn.metrics.cluster import adjusted_rand_score, silhouette_score
import collections
from sklearn.cluster import KMeans

base_path = '/Users/zhongyuanke/data/'
vcca_path = base_path + 'vcca/pbmc_protein/vcca_temp_04.h5ad'
seurat_path = base_path + 'vcca/pbmc_protein/seurat_protein.h5ad'
totalvi_path = base_path + 'vcca/pbmc_protein/totalvi/protein_totalvi.h5ad'


adata_vcca = sc.read_h5ad(vcca_path)
adata_seurat = sc.read_h5ad(seurat_path)
adata_totalvi = sc.read_h5ad(totalvi_path)

celltype = adata_vcca.obs['vcca_celltype_l2']

print(collections.Counter(celltype))
label = tl.text_label_to_number(celltype)
print(collections.Counter(label))


my_color_map = {
    0: 'deeppink',
    1: 'royalblue',
    2: 'teal',
    3: 'limegreen',
    4: 'grey',
    5: 'olive',
    6: 'skyblue',
    7: 'darkslateblue',
    8: 'violet',
    9: 'blueviolet',
    10: 'sandybrown',
    11: 'gold',
    12: 'lightcoral',
    13: 'turquoise',
    14: 'slateblue',
    15: 'rosybrown',
    16: 'red',
    17: 'darkred',
    18: 'tan',
    19: 'saddlebrown',
    20: 'deepskyblue',
    21: 'cyan',
    22: 'green',
    23: 'brown',
    24: 'orange',
    25: 'fuchsia'
}
label_color = []
for l in label:
    label_color.append(my_color_map[l])


print(adata_seurat)
print(adata_totalvi)
data_seurat_emb = adata_seurat.obsm['X_wnn.umap']

sc.pp.neighbors(adata_totalvi,use_rep='X_totalVI')
sc.tl.umap(adata_totalvi)
data_totalvi_emb = adata_totalvi.obsm['X_umap']

data_vimcca_emb = adata_vcca.obsm['X_umap']

print(data_totalvi_emb.shape)

batch_cmap = 'Dark2'
c_map = 'tab20'
size = 1
xy_label_size = 10
title_size = 10
xy_label = 'UMAP_'
alpha=0.8

fig = plt.figure(figsize=(50, 30))
plt.subplots_adjust(wspace =0.3, hspace =0.3)
ax = fig.add_subplot(131)
plt.title('Seurat 4.0', fontsize=title_size)
ax.scatter(data_seurat_emb[:, 0], data_seurat_emb[:, 1], c=label_color, cmap=batch_cmap, s=size,
           linewidth=0,alpha=alpha)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel(xy_label + '1', size=xy_label_size)
plt.ylabel(xy_label + '2', size=xy_label_size)


ax = fig.add_subplot(132)
plt.title('TotalVI', fontsize=title_size)
ax.scatter(data_totalvi_emb[:, 0], data_totalvi_emb[:, 1], c=label_color, cmap=c_map, s=size,
           linewidth=0,alpha=alpha)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel(xy_label + '1', size=xy_label_size)
plt.ylabel(xy_label + '2', size=xy_label_size)

ax = fig.add_subplot(133)
plt.title('VIMCCA', fontsize=title_size)
ax.scatter(data_vimcca_emb[:, 0], data_vimcca_emb[:, 1], c=label_color, cmap=c_map, s=size,
           linewidth=0,alpha=alpha)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel(xy_label + '1', size=xy_label_size)
plt.ylabel(xy_label + '2', size=xy_label_size)

plt.show()
plt.close(fig)
# --------------------ARI-------------------------------
print(len(set(label)))
cluster_count = 21
kmeans_seurat = KMeans(n_clusters=cluster_count).fit(adata_seurat.obsm['X_wnn.umap'])
ari_seurat = adjusted_rand_score(label, kmeans_seurat.labels_)
sh_seurat = adjusted_rand_score(label, kmeans_seurat.labels_)
print('Seurat 4: ', ari_seurat, sh_seurat)

kmeans_totalvi = KMeans(n_clusters=cluster_count).fit(adata_totalvi.obsm['X_totalVI'])
ari_totalvi = adjusted_rand_score(label, kmeans_totalvi.labels_)
sh_totalvi = adjusted_rand_score(label, kmeans_totalvi.labels_)
print('TotalVI: ', ari_totalvi, sh_totalvi)

kmeans_vcca = KMeans(n_clusters=cluster_count).fit(adata_vcca.obsm['X_vcca'])
ari_vcca = adjusted_rand_score(label, kmeans_vcca.labels_)
sh_vcca = adjusted_rand_score(label, kmeans_vcca.labels_)
print('VCCA: ', ari_vcca, sh_vcca)



