import sys
import argparse
import scvi
import scanpy as sc

parser = argparse.ArgumentParser()
parser.add_argument("--base_path", type=str, default='/Users/zhongyuanke/data/', help="base path")
opt = parser.parse_args()


base_path = opt.base_path
file_rna = base_path + 'multimodal/protein/vdj_v1_hs_pbmc2_5gex_protein_filtered_feature_bc_matrix.h5'
file_protein = base_path + 'multimodal/protein/protein.csv'

adata= sc.read_10x_h5(file_rna)
adata_y = sc.read_csv(file_protein)
sc.pp.filter_genes(adata, min_cells=3)
#
print(adata)
adata.obsm['protein']=adata_y.X
adata.layers["counts"]=adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
scvi.model.TOTALVI.setup_anndata(
    adata,
    protein_expression_obsm_key="protein",
    layer="counts",
)

# adata = scvi.data.pbmcs_10x_cite_seq(run_setup_anndata=False)
# print(adata)
# print(adata.obs['batch'])
# adata.layers["counts"]=adata.X.copy()
# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata)
# adata.raw = adata
# #
# scvi.model.TOTALVI.setup_anndata(
#     adata,
#     protein_expression_obsm_key="protein_expression",
#     layer="counts",
#     batch_key="batch"
# )

print(adata)
vae = scvi.model.TOTALVI(adata, latent_distribution="normal")
vae.train(max_epochs=200)
adata.obsm["X_totalVI"] = vae.get_latent_representation()
adata.write_h5ad('/Users/zhongyuanke/data/vcca/pbmc_protein/totalvi/protein_totalvi.h5ad')