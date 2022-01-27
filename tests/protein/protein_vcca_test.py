import scanpy as sc
from scbean.model import vcca as vcca
from scbean.tools import utils as tl
import matplotlib
import pandas as pd
matplotlib.use('TkAgg')

base_path = '/Users/zhongyuanke/data/'
file_rna = base_path + 'multimodal/protein/rna.csv'
file_rna = base_path + 'multimodal/protein/vdj_v1_hs_pbmc2_5gex_protein_filtered_feature_bc_matrix.h5'
file_protein = base_path + 'multimodal/protein/protein.csv'

adata_x= sc.read_10x_h5(file_rna)
adata_y = sc.read_csv(file_protein)


sc.pp.filter_genes(adata_x, min_cells=1)
sc.pp.log1p(adata_x)
sc.pp.log1p(adata_y)
sc.pp.scale(adata_x)
sc.pp.scale(adata_y)
# sc.pp.highly_variable_genes(adata_x, n_top_genes=5000)
# adata_x = adata_x[:, adata_x.var.highly_variable]

z = vcca.fit_integration(
    adata_x, adata_y,
    sparse_x=False,
    sparse_y=False,
    mode='VCCA',
    hidden_layers=[128, 64, 32, 8],
    epochs=50
)

adata_x.obsm['X_vcca'] = z
sc.pp.neighbors(adata_x, use_rep='X_vcca', n_neighbors=13)
sc.tl.leiden(adata_x, resolution=2)
sc.tl.louvain(adata_x, resolution=2)
sc.tl.umap(adata_x)
# sc.pl.umap(adata_x, color=['CST3', 'CD8A', 'CD79A', 'CD74', 'MS4A1', 'LTB', 'FCER1G', 'TYROBP',
#                            'CD3E', 'CCL5', 'TYROBP', 'TCF7', 'LYZ', 'IL7R', 'HLA-DRA', 'LDHB',  'louvain'], s=2, ncols=4)
sc.pl.umap(adata_x, color='leiden', legend_loc='on data')

adata_x.write_h5ad('/Users/zhongyuanke/data/vcca/pbmc_protein/vcca_temp_04.h5ad')

dataframe = pd.DataFrame(adata_x.obs['leiden'])
dataframe.to_csv("/Users/zhongyuanke/data/vcca/pbmc_protein/vcca_leiden.csv",index=False,sep=',')
print('finishi write leiden')

dataframe = pd.DataFrame(adata_x.obs['louvain'])
dataframe.to_csv("/Users/zhongyuanke/data/vcca/pbmc_protein/vcca_louvain.csv",index=False,sep=',')
print('finishi write louvain')
