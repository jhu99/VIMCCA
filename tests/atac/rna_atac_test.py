import scanpy as sc
from network.vcca import PVCCA
from scbean.model import vcca
from utils import tools as tl
import pandas as pd
import umap
import numpy as np
import matplotlib
matplotlib.use('TkAgg')


base_path = '/Users/zhongyuanke/data/'
file_rna = base_path + 'multimodal/atac_pbmc_10k/rna_filt.csv'
file_atac = base_path + 'multimodal/atac_pbmc_10k/atac.h5ad'
file_activity = base_path + 'multimodal/atac_pbmc_10k/activaty.h5ad'
celltype_path = base_path + 'multimodal/atac_pbmc_10k/seurat_celltype_filt.csv'
batch_size = 128

adata_x = sc.read_csv(file_rna)

adata_y = sc.read_h5ad(file_atac)

sc.pp.filter_genes(adata_x, min_cells=1)
sc.pp.filter_genes(adata_y, min_cells=1)
sc.pp.log1p(adata_x)
sc.pp.log1p(adata_y)
sc.pp.scale(adata_x)
sc.pp.scale(adata_y)

# adata_x.write_h5ad('/Users/zhongyuanke/data/vcca/atac/pbmc_10k/pvcca_temp.h5ad')
celltype = tl.get_label_by_txt(celltype_path)
label = tl.text_label_to_number(celltype)
adata_x.obs['seurat_celltype'] = celltype
# adata_x.write_h5ad('/Users/zhongyuanke/data/multimodal/atac_pbmc_10k/rna.h5ad')
print(adata_y)


# sc.pp.scale(adata_x)

# adata_x.obs['label'] = label
# adata_x.obs['celltype'] = celltype
# sc.pp.neighbors(adata_x)
# sc.tl.louvain(adata_x)
# sc.tl.umap(adata_x)
# sc.pl.umap(adata_x, color='celltype', color_map='tab20c',sort_order=False, legend_loc='right margin', s=2)
#
# adata_y.obs['label'] = label
# adata_y.obs['celltype'] = celltype
# sc.pp.neighbors(adata_y)
# sc.tl.louvain(adata_y)
# sc.tl.umap(adata_y)
# sc.pl.umap(adata_y, color='celltype', color_map='tab20c', sort_order=False, legend_loc='right margin', s=2)

# sc.pp.log1p(adata_y)
# sc.pp.highly_variable_genes(adata_x, n_top_genes=5000)
# sc.pp.highly_variable_genes(adata_y, n_top_genes=5000)
print(adata_x)
print(adata_y)


z = vcca.fit_integration(
    adata_x,
    adata_y,
    sparse_x=False,
    sparse_y=False,
    mode='VCCA',
    hidden_layers=[128, 64, 16, 8],
    epochs=50)
# adata_x.write_h5ad(base_path+'vcca/atac/pbmc_10k/pvcca_temp.h5ad')
# z=np.array(z)
# print(z.shape[1])
# z_emb = umap.UMAP().fit_transform(z)
adata_x.obsm['X_vcca'] = z
sc.pp.neighbors(adata_x, use_rep='X_vcca', n_neighbors=13)
# sc.tl.louvain(adata_x)
sc.tl.leiden(adata_x, resolution=2)
sc.tl.umap(adata_x)
adata_x.obs['seurat_celltype'] = celltype
sc.pl.umap(adata_x, color='seurat_celltype', legend_loc='right margin', color_map='tab20c', s=2)
print(adata_x)

# adata_x.obsm['X_umap'] = z_emb
# adata_x.obs['celltype'] = celltype
adata_x.write_h5ad(base_path+'vcca/atac/pbmc_10k/vcca_temp.h5ad')


dataframe = pd.DataFrame(adata_x.obs['leiden'])
dataframe.to_csv('/Users/zhongyuanke/data/vcca/atac/pbmc_10k/vcca_leiden.csv',index=False,sep=',')
print('finishi write leiden')

# sc.pl.umap(adata_x, color=['CD79A', 'S100A8', 'GNLY', 'CST3', 'CD3D', 'celltype'], legend_loc='right margin')
# sc.pl.umap(adata_x, color='celltype', legend_loc='right margin', color_map='tab20c', s=2)
# print(adata_x)
# pred_label = adata_x.obs['louvain'].values
# pred_label = np.array(pred_label, dtype=int)

# z_emb = umap.UMAP().fit_transform(z)
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.scatter(z_emb[:, 0], z_emb[:, 1], c=label, cmap='tab20b', s=1, linewidth=0)
# plt.show()


