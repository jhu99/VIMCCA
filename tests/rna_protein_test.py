import scanpy as sc
# from network.vcca import VCCA
from network.p_vcca import VCCA
from utils import tools as tl
import umap
import numpy as np
import matplotlib.pyplot as plt


base_path = '/Users/zhongyuanke/data/'
file_rna = base_path + 'multimodal/bmcite/rna_5k.csv'
file_atac = base_path + 'multimodal/bmcite/protein.csv'
celltype_path = base_path + 'multimodal/bmcite/celltype_l2.csv'

adata_x = sc.read_csv(file_rna)
adata_y = sc.read_csv(file_atac)

celltype = tl.get_label_by_txt(celltype_path)
label = tl.text_label_to_number(celltype)
print(label)
adata_x.obs['celltype_l1'] = celltype
adata_x.obs['label'] = label
#
# sc.pl.umap(adata_x, color=['CST3', 'CD8A', 'CD79A', 'CD74',
#                            'CCL5', 'TYROBP', 'LYZ', 'louvain', 'celltype_l1'],
#            s=2)
adata_y.obs['celltype_l1'] = celltype

# sc.pp.neighbors(adata_y)
# sc.tl.louvain(adata_y)
# sc.tl.umap(adata_y)
# sc.pl.umap(adata_y, color=["CD123", "CD127-IL7Ra", "CD14", "CD16", 'louvain', 'celltype_l1'])
# sc.pl.umap(adata_y, color='celltype_l1', legend_loc='on data')

# sc.pp.filter_genes(adata_x, min_cells=30)
# sc.pp.filter_genes(adata_y, min_cells=30)
# sc.pp.log1p(adata_x)
# sc.pp.log1p(adata_y)
# sc.pp.highly_variable_genes(adata_x, n_top_genes=5000)
# sc.pp.highly_variable_genes(adata_y, n_top_genes=5000)
print(adata_x)
print(adata_y)
# adata_x = adata_x[:adata_x.var['highly_variable']]

x = adata_x.X
y = adata_y.X

net_x = VCCA(input_size_x=x.shape[1], inputs_size_y=y.shape[1], latent_z_size=6, latent_xy_size=2)
net_x.build()
net_x.compile()
his = net_x.train(x, y, epochs=45, batch_size=batch_size)
z = net_x.integrate(x)

# x_emb = umap.UMAP().fit_transform(x)

z_emb = umap.UMAP().fit_transform(z)
adata_x.obsm['vcca'] = z
adata_x.obs['label'] = label


adata_y.obsm['vcca'] = z
adata_y.obs['label'] = label
adata_y.obs['celltype_l1'] = celltype

sc.pp.neighbors(adata_x, use_rep='vcca')
# sc.pp.neighbors(adata_x)
sc.tl.louvain(adata_x)
sc.tl.umap(adata_x)
adata_x.obsm['X_umap'] = z_emb
adata_x.write_h5ad(base_path+'vcca/bmcite/pvcca_01.h5ad')
# sc.pl.umap(adata_x, color=['CST3', 'CD8A', 'CD79A', 'CD74',
#                            'CCL5', 'TYROBP', 'LYZ',  'louvain', 'celltype_l1'],
#            s=2)
sc.pl.umap(adata_x, color='celltype_l1', legend_loc='on data')


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


