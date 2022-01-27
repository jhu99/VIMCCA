from sklearn.metrics.cluster import adjusted_rand_score, silhouette_score
from sklearn.cluster import KMeans
import scanpy as sc
import pandas as pd


base_path = '/Users/zhongyuanke/data/'
wnn_umap_path='/Users/zhongyuanke/data/multimodal/atac_pbmc_10k/wnn_umap_filt.csv'
wnn_umap = pd.read_csv(wnn_umap_path, index_col=0)
wnn_umap = wnn_umap.values
file_rna = '/Users/zhongyuanke/data/multimodal/atac_pbmc_10k/rna.h5ad'
file_vcca = base_path + 'vcca/atac/pbmc_10k/vcca_temp.h5ad'
adata_vcca = sc.read_h5ad(file_vcca)

vcca_label = adata_vcca.obs['vcca_celltype_l2']
kmeans_seurat = KMeans(n_clusters=15).fit(wnn_umap)
ari_seurat = adjusted_rand_score(vcca_label, kmeans_seurat.labels_)
print('ARI_seurat  ',ari_seurat)

kmeans_vcca = KMeans(n_clusters=15).fit(adata_vcca.obsm['X_umap'])
ari_vcca = adjusted_rand_score(vcca_label, kmeans_vcca.labels_)
print('ARI_VCCA  ',ari_vcca)