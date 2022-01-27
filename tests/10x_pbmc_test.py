import scanpy as sc
from network.vcca import VCCA, PVCCA
from utils import tools as tl
import umap
import numpy as np
import matplotlib.pyplot as plt


base_path = '/Users/zhongyuanke/data/'

file_rna = base_path + 'multimodal/protein/rna.csv'
file_atac = base_path + 'multimodal/protein/protein.csv'

batch_size = 128

adata_x = sc.read_csv(file_rna)
adata_y = sc.read_csv(file_atac)




sc.pp.filter_genes(adata_x, min_cells=80)
# sc.pp.filter_genes(adata_y, min_cells=30)
sc.pp.log1p(adata_x)
sc.pp.log1p(adata_y)
sc.pp.scale(adata_x)
sc.pp.scale(adata_y)
# sc.pp.highly_variable_genes(adata_x, n_top_genes=8000)
# sc.pp.highly_variable_genes(adata_y, n_top_genes=8000)
# adata_x = adata_x[:, adata_x.var.highly_variable]
# adata_y = adata_y[:, adata_y.var.highly_variable]
print(adata_x)
print(adata_y)
# sc.pp.neighbors(adata_x)
# sc.tl.louvain(adata_x)
# sc.tl.umap(adata_x)
# sc.pl.umap(adata_x, color=['S100A8', 'CST3', 'CD8A', 'CD79A', 'CD74',
#                            'CD3E', 'CCL5', 'TYROBP', 'TCF7', 'LYZ', 'IL7R', 'louvain'], s=2, ncols=4)
#
# sc.pp.neighbors(adata_y)
# sc.tl.louvain(adata_y)
# sc.tl.umap(adata_y)
# sc.pl.umap(adata_y, color=['CD3_TotalSeqC', 'CD45RA_TotalSeqC', 'CD8a_TotalSeqC', 'CD16_TotalSeqC',
#                            'CD25_TotalSeqC', 'CD45RO_TotalSeqC', 'PD-1_TotalSeqC', 'TIGIT_TotalSeqC',
#                            'CD127_TotalSeqC', 'CD15_TotalSeqC', 'louvain'], s=2, ncols=4)
x = adata_x.X
y = adata_y.X

net_x = PVCCA(input_size_x=x.shape[1], inputs_size_y=y.shape[1], hidden_layers=[128, 64, 32, 5])
net_x.build()
net_x.compile()
his = net_x.train(x, y, epochs=30, batch_size=batch_size)
z = net_x.integrate(x)
adata_x.obsm['vcca'] = z

# x_emb = umap.UMAP().fit_transform(x)
# z_emb = umap.UMAP().fit_transform(z)

sc.pp.neighbors(adata_x, use_rep='vcca')
# sc.pp.neighbors(adata_x)
sc.tl.louvain(adata_x)
sc.tl.umap(adata_x)
# adata_x.write_h5ad(base_path+'vcca/pbmc_protein/pvcca_03.h5ad')
sc.tl.rank_genes_groups(adata_x, 'louvain', method='t-test')
sc.pl.rank_genes_groups(adata_x, n_genes=5, sharey=False)
# adata_x.obsm['X_umap'] = z_emb
sc.pl.umap(adata_x, color=['S100A8', 'CST3', 'CD8A', 'CD79A', 'CD74', 'MS4A1', 'LTB', 'FCER1G', 'TYROBP',
                           'CD3E', 'CCL5', 'TYROBP', 'TCF7', 'LYZ', 'IL7R', 'HLA-DRA', 'LDHB',  'louvain'], s=2, ncols=4)
sc.pl.umap(adata_x, color='louvain', legend_loc='on data')
# sc.pl.umap(adata_x, color='louvain', legend_loc='on data')
# sc.pl.umap(adata_y, color=["CD123", "CD127-IL7Ra", "CD14", "CD16", 'CD79A', 'S100A8', 'HBB', 'GNLY', 'CD3D', 'louvain'])

# sc.pp.neighbors(adata_y, use_rep='vcca')
# sc.pp.neighbors(adata_y)
# sc.tl.louvain(adata_y)
# sc.tl.umap(adata_y)
# sc.pl.umap(adata_y, color=["CD123", "CD127-IL7Ra", "CD14", "CD16", 'louvain', 'celltype_l1'])

# pred_label = adata_x.obs['louvain'].values
# pred_label = np.array(pred_label, dtype=int)

# z_emb = umap.UMAP().fit_transform(z)
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.scatter(z_emb[:, 0], z_emb[:, 1], c=label, cmap='tab20b', s=1, linewidth=0)
# plt.show()


# fig = plt.figure()
# ax = fig.add_subplot(121)
# ax.scatter(x_emb[:, 0], x_emb[:, 1], c=label, cmap='tab20b', s=2, linewidth=0)
#
# y_emb = umap.UMAP().fit_transform(y)
# ax = fig.add_subplot(122)
# ax.scatter(y_emb[:, 0], y_emb[:, 1], c=label, cmap='tab20b', s=1, linewidth=0)
# plt.show()


