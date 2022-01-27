import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from utils import tools as tl
import utils.tools as tt
from sklearn.preprocessing import LabelEncoder

base_path = '/Users/zhongyuanke/data/'

file_rna = base_path + 'multimodal/protein/vdj_v1_hs_pbmc2_5gex_protein_filtered_feature_bc_matrix.h5'
file_protein = base_path + 'multimodal/protein/protein.csv'
# file_vcca = base_path + 'vcca/pbmc_protein/vcca_temp.h5ad'

# adata_vcca = sc.read_h5ad(file_vcca)
rna = sc.read_10x_h5(file_rna)
protein = sc.read_csv(file_protein)

# rna.obs['leiden']=adata_vcca.obs['leiden']
# rna.obs['vcca_celltype']=adata_vcca.obs['vcca_celltype']

# sc.pp.log1p(rna)
# sc.pp.pca(rna)
# sc.pp.neighbors(rna)
# sc.tl.umap(rna)
#
#
# sc.pp.log1p(protein)
# sc.pp.pca(protein)
# sc.pp.neighbors(protein)
# sc.tl.umap(protein)


# protein.obs['leiden']=adata_vcca.obs['leiden']
# protein.obs['vcca_celltype']=adata_vcca.obs['vcca_celltype']
rna.obs_names_make_unique()
rna.var_names_make_unique()

print(rna)
rna.write_h5ad('/Users/zhongyuanke/data/vcca/pbmc_protein/R/rna_orig.h5ad')
protein.write_h5ad('/Users/zhongyuanke/data/vcca/pbmc_protein/R/protein_orig.h5ad')
# adata_vcca.write_h5ad('/Users/zhongyuanke/data/vcca/pbmc_protein/R/vcca.h5ad')