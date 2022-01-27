import scanpy as sc
import pandas as pd
from matplotlib.pyplot import rc_context

base_path = '/Users/zhongyuanke/data/'
file_vcca = base_path + 'vcca/atac/pbmc_10k/vcca_temp.h5ad'

adata = sc.read_h5ad(file_vcca)

keys1 = ['CAMK4',] # CAMK4(naive CD4 T-cell, MAIT T-cell, T-reg, naive CD8 T-cell, memory CD4 T-cell, memory CD8 T-cell)
keys2 = ['CD79A'] # naive B, memory B
keys3 = ['ZFAT', 'TCF4'] # ZFAT(plasmacytoid DC) TCF4 (plasmacytoid DC, memory B-cell, naive B-cell)
keys4 = ['FCER1G','LYPD2'] # FCER1G(non-classical monocyte)
keys5 = ['LEF1',] # BCL11b(naive CD4 T-cell, naive CD8 T-cell, memory CD4 T-cell, memory CD8 T-cell, gdT-cell, T-reg, MAIT T-cell)
keys6 = ['STAM']

with rc_context({'figure.figsize': (6.5, 3)}):
    sc.pl.violin(adata, keys=keys4, groupby='vcca_celltype_l2', stripplot=False, inner='box')
    sc.pl.violin(adata, keys=keys4, groupby='seurat_celltype', stripplot=False, inner='box')


# sc.pl.(adata, ['CAMK4', 'FCER1G'],groupby='vcca_celltype_l2',
#              jitter=0.4, multi_panel=True)
# gene_names = pd.Series(adata.var_names, index=adata.var_names)
# lookup_keys = []
# not_found = []
# keys = ['CAMK4', 'FCER1G','LYPD2','STAM','CST7']
# for key in keys:
#     if key in adata.obs.columns:
#         lookup_keys.append(key)
#     elif key in gene_names.index:
#         lookup_keys.append(gene_names[key])
#     else:
#         not_found.append(key)
# df = pd.DataFrame(index=adata.obs_names)
# for k, l in zip(keys, lookup_keys):
#     df[k] = adata.raw.obs_vector(l)
#
# print(df)
