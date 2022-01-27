import scanpy as sc
from scbean.model import vcca as vcca
from scbean.tools import utils as tl
import matplotlib
import pandas as pd
import time
matplotlib.use('TkAgg')

base_path = '/Users/zhongyuanke/data/'
file_rna = base_path + 'multimodal/protein/rna.csv'
file_rna = base_path + 'multimodal/protein/vdj_v1_hs_pbmc2_5gex_protein_filtered_feature_bc_matrix.h5'
file_protein = base_path + 'multimodal/protein/protein.csv'

adata_x= sc.read_10x_h5(file_rna)
adata_y = sc.read_csv(file_protein)

adata2=adata_x
adata2.var_names_make_unique(join='2')
adata_x = adata_x.concatenate(adata2)

adata_y2 = adata_y
adata_y2.var_names_make_unique(join='2')
adata_y = adata_y2.concatenate(adata_y2)


# --------------------3000---------------------------
# adata_x_temp = adata_x[0:3000]
# adata_y_temp = adata_y[0:3000]
# t0 = time.time()
#
# sc.pp.filter_genes(adata_x_temp, min_cells=1)
# sc.pp.log1p(adata_x_temp)
# sc.pp.log1p(adata_y_temp)
# sc.pp.scale(adata_x_temp)
# sc.pp.scale(adata_y_temp)
#
# z = vcca.fit_integration(
#     adata_x_temp, adata_y_temp,
#     sparse_x=False,
#     sparse_y=False,
#     mode='VCCA',
#     hidden_layers=[128, 64, 32, 8],
#     epochs=35
# )
#
# t1 = time.time()
# print("Total time running VCCA 3000 cells : %s min" % (str((t1-t0)/60)))
#
# # --------------------6000---------------------------
# adata_x_temp = adata_x[0:6000]
# adata_y_temp = adata_y[0:6000]
# t0 = time.time()
#
# sc.pp.filter_genes(adata_x_temp, min_cells=1)
# sc.pp.log1p(adata_x_temp)
# sc.pp.log1p(adata_y_temp)
# sc.pp.scale(adata_x_temp)
# sc.pp.scale(adata_y_temp)
#
# z = vcca.fit_integration(
#     adata_x_temp, adata_y_temp,
#     sparse_x=False,
#     sparse_y=False,
#     mode='VCCA',
#     hidden_layers=[128, 64, 32, 8],
#     epochs=35
# )
#
# t1 = time.time()
# print("Total time running VCCA 6000 cells : %s min" % (str((t1-t0)/60)))
#
# # --------------------------9000----------------------------------
# adata_x_temp = adata_x[0:9000]
# adata_y_temp = adata_y[0:9000]
# t0 = time.time()
#
# sc.pp.filter_genes(adata_x_temp, min_cells=1)
# sc.pp.log1p(adata_x_temp)
# sc.pp.log1p(adata_y_temp)
# sc.pp.scale(adata_x_temp)
# sc.pp.scale(adata_y_temp)
#
# z = vcca.fit_integration(
#     adata_x_temp, adata_y_temp,
#     sparse_x=False,
#     sparse_y=False,
#     mode='VCCA',
#     hidden_layers=[128, 64, 32, 8],
#     epochs=35
# )
#
# t1 = time.time()
# print("Total time running VCCA 9000 cells : %s min" % (str((t1-t0)/60)))
#
# # --------------------------12000----------------------------------
adata_x_temp = adata_x[0:12000]
adata_y_temp = adata_y[0:12000]
t0 = time.time()

sc.pp.filter_genes(adata_x_temp, min_cells=1)
sc.pp.log1p(adata_x_temp)
sc.pp.log1p(adata_y_temp)
sc.pp.scale(adata_x_temp)
sc.pp.scale(adata_y_temp)

z = vcca.fit_integration(
    adata_x_temp, adata_y_temp,
    sparse_x=False,
    sparse_y=False,
    mode='VCCA',
    hidden_layers=[128, 64, 32, 8],
    epochs=35
)

t1 = time.time()
print("Total time running VCCA 12000 cells : %s min" % (str((t1-t0)/60)))

# --------------------------15000----------------------------------
# adata_x_temp = adata_x[0:15000]
# adata_y_temp = adata_y[0:15000]
# t0 = time.time()
#
# sc.pp.filter_genes(adata_x_temp, min_cells=1)
# sc.pp.log1p(adata_x_temp)
# sc.pp.log1p(adata_y_temp)
# sc.pp.scale(adata_x_temp)
# sc.pp.scale(adata_y_temp)
#
# z = vcca.fit_integration(
#     adata_x_temp, adata_y_temp,
#     sparse_x=False,
#     sparse_y=False,
#     mode='VCCA',
#     hidden_layers=[128, 64, 32, 8],
#     epochs=30
# )
#
# t1 = time.time()
# print("Total time running VCCA 15000 cells : %s min" % (str((t1-t0)/60)))