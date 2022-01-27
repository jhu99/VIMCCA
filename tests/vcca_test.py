import scanpy as sc
import network.vcca as vcca
from utils import tools as tl
import matplotlib
matplotlib.use('TkAgg')

base_path = '/Users/zhongyuanke/data/'
adata_path = base_path + 'multimodal/protein/vdj_v1_hs_pbmc2_5gex_protein_filtered_feature_bc_matrix.h5'
file_rna = base_path + 'multimodal/pbmc_10k/rna_norm.csv'
file_atac = base_path + 'multimodal/pbmc_10k/lsi.csv'
celltype_path = base_path + 'multimodal/pbmc_10k/celltype_filt.csv'
batch_size = 128
adata = sc.read_h5ad(adata_path)
adata_x = sc.read_csv(file_rna)
adata_y = sc.read_csv(file_atac)
print(adata)
celltype = tl.get_label_by_txt(celltype_path)
label = tl.text_label_to_number(celltype)
print(label)

sc.pp.filter_genes(adata_x, min_cells=100)
sc.pp.filter_genes(adata_y, min_cells=30)
print(adata_x)
print(adata_y)


z = vcca.fit_integration(adata_x, adata_y, hidden_layers=[128, 64,], latent_xy_size=2)

adata_x.obsm['X_vcca'] = z
adata_x.obs['celltype'] = celltype

sc.pp.neighbors(adata_x, use_rep='X_vcca')
sc.tl.louvain(adata_x)
sc.tl.umap(adata_x)



