import sys
import argparse
import scvi
import scanpy as sc
import time


parser = argparse.ArgumentParser()
parser.add_argument("--base_path", type=str, default='/Users/zhongyuanke/data/', help="base path")
opt = parser.parse_args()


base_path = opt.base_path
file_rna = base_path + 'multimodal/protein/vdj_v1_hs_pbmc2_5gex_protein_filtered_feature_bc_matrix.h5'
file_protein = base_path + 'multimodal/protein/protein.csv'

adata_x= sc.read_10x_h5(file_rna)
adata_y = sc.read_csv(file_protein)
sc.pp.filter_genes(adata_x, min_cells=1)

adata2=adata_x
adata2.var_names_make_unique(join='2')
adata_x = adata_x.concatenate(adata2)
adata_x.var_names_make_unique()

adata_y2 = adata_y
adata_y2.var_names_make_unique(join='2')
adata_y = adata_y2.concatenate(adata_y2)
adata_y.var_names_make_unique()


# ----------------------------3000------------------------------------
adata_x_temp= adata_x[0:3000]
adata_y_temp = adata_y[0:3000]
print(adata_x_temp)


t0=time.time()

adata_x_temp.obsm['protein']=adata_y_temp.X
adata_x_temp.layers["counts"]=adata_x_temp.X.A
sc.pp.normalize_total(adata_x_temp, target_sum=1e4)
sc.pp.log1p(adata_x_temp)
adata_x_temp.raw = adata_x_temp
scvi.model.TOTALVI.setup_anndata(
    adata_x_temp,
    protein_expression_obsm_key="protein",
    layer="counts",
)
# print(dict(adata_x_temp.attrs))
print(adata_x_temp)
vae = scvi.model.TOTALVI(adata_x_temp, latent_distribution="normal")
vae.train(max_epochs=200)
adata_x_temp.obsm["X_totalVI"] = vae.get_latent_representation()
t1 = time.time()
print("Total time running TotalVI on 3000 cell : %s seconds" % (str(t1-t0)))


# ----------------------------6000------------------------------------
# adata_1 = adata[0:6000]
# adata_y_1 = adata_y[0:6000]
#
# t0=time.time()
# sc.pp.filter_genes(adata_1, min_cells=1)
# adata.obsm['protein']=adata_y_1.X
# adata_1.layers["counts"]=adata_1.X.copy()
# sc.pp.normalize_total(adata_1, target_sum=1e4)
# sc.pp.log1p(adata_1)
# adata_1.raw = adata_1
# scvi.model.TOTALVI.setup_anndata(
#     adata_1,
#     protein_expression_obsm_key="protein",
#     layer="counts",
# )
#
# vae = scvi.model.TOTALVI(adata_1, latent_distribution="normal")
# vae.train(max_epochs=200)
# adata_1.obsm["X_totalVI"] = vae.get_latent_representation()
# t1 = time.time()
# print("Total time running TotalVI on 6000 cell : %s seconds" % (str(t1-t0)))

# ----------------------------12000------------------------------------
# adata_1 = adata[0:12000]
# adata_y_1 = adata_y[0:12000]
#
# t0 = time.time()
# sc.pp.filter_genes(adata_1, min_cells=1)
# adata.obsm['protein']=adata_y_1.X
# adata_1.layers["counts"]=adata_1.X.copy()
# sc.pp.normalize_total(adata_1, target_sum=1e4)
# sc.pp.log1p(adata_1)
# adata_1.raw = adata_1
# scvi.model.TOTALVI.setup_anndata(
#     adata_1,
#     protein_expression_obsm_key="protein",
#     layer="counts",
# )
# vae = scvi.model.TOTALVI(adata_1, latent_distribution="normal")
# vae.train(max_epochs=200)
# adata_1.obsm["X_totalVI"] = vae.get_latent_representation()
# t1 = time.time()
# print("Total time running TotalVI on 12000 cell : %s seconds" % (str(t1-t0)))