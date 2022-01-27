import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from utils import tools as tl
import utils.tools as tt
from sklearn.preprocessing import LabelEncoder

base_path = '/Users/zhongyuanke/data/'

file_rna = base_path + 'multimodal/protein/vdj_v1_hs_pbmc2_5gex_protein_filtered_feature_bc_matrix.h5'
file_protein = base_path + 'multimodal/protein/protein.csv'
file_vcca = base_path + 'vcca/pbmc_protein/vcca_temp_04.h5ad'

adata_vcca = sc.read_h5ad(file_vcca)

rna = sc.read_10x_h5(file_rna)
# sc.pp.filter_genes(rna, min_cells=30)
sc.pp.log1p(rna)
sc.pp.pca(rna)
sc.pp.neighbors(rna)
sc.tl.umap(rna)
rna.obs['vcca_celltype'] = adata_vcca.obs['vcca_celltype_l2']
# sc.pl.umap(rna, color='vcca_celltype', legend_loc='on data')
# sc.pl.umap(rna, color='vcca_leiden')
# sc.pl.umap(rna, color = 'seurat_celltype')

# ---------------------------------------------------------------------------
protein = sc.read_csv(file_protein)
sc.pp.log1p(protein)
# sc.pp.filter_genes(atac, min_cells=30)
sc.pp.pca(protein)
sc.pp.neighbors(protein)
sc.tl.umap(protein)


rna_emb = rna.obsm['X_umap']
atac_emb = protein.obsm['X_umap']
vcca_emb = adata_vcca.obsm['X_umap']

celltype = adata_vcca.obs['vcca_celltype_l2'].values
encoder = LabelEncoder()

label = encoder.fit_transform(celltype)
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

# label = list(map(int, adata_vcca.obs['leiden']))
# seurat_celltype_path = base_path + 'multimodal/atac_pbmc_10k/seurat_celltype_filt.csv'
#
# label = tl.get_label_by_txt(seurat_celltype_path)
# label = encoder.fit_transform(celltype)
# c_map = 'Dark2'
c_map = 'tab20c'
size = 2
xy_label_size = 10
title_size = 14
xy_label = 'UMAP_'

fig = plt.figure(figsize=(50, 30))
plt.subplots_adjust(wspace =0.3, hspace =0.3)

ax = fig.add_subplot(131)
plt.title('RNA', fontsize=title_size)
ax.scatter(rna_emb[:, 0], rna_emb[:, 1], c=label_color, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel(xy_label + '1', size=xy_label_size)
plt.ylabel(xy_label + '2', size=xy_label_size)


ax = fig.add_subplot(132)
plt.title('Protein', fontsize=title_size)
ax.scatter(atac_emb[:, 0], atac_emb[:, 1], c=label_color, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel(xy_label + '1', size=xy_label_size)
plt.ylabel(xy_label + '2', size=xy_label_size)

ax = fig.add_subplot(133)
plt.title('VCCA', fontsize=title_size)
ax.scatter(vcca_emb[:, 0], vcca_emb[:, 1], c=label_color, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel(xy_label + '1', size=xy_label_size)
plt.ylabel(xy_label + '2', size=xy_label_size)
plt.show()
plt.close(fig)