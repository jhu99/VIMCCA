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
file_vcca = base_path + 'vcca/atac/pbmc_10k/vcca_temp.h5ad'
celltype_path = base_path + 'multimodal/atac_pbmc_10k/seurat_celltype_filt.csv'

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
adata.obs['seurat_celltype']=label
# adata.write_h5ad(file_vcca)

sc.pp.log1p(adata)
sc.tl.louvain(adata, resolution=1.6)

adata0 = adata[adata.obs['louvain']=='0']

sc.pp.log1p(adata0)
sc.tl.louvain(adata0, resolution=1.6)
sc.tl.umap(adata0)
sc.pl.umap(adata0, color='louvain',legend_loc='on data')
print(adata0)