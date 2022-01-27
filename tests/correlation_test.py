import scanpy as sc
import numpy as np
import anndata
from sklearn.preprocessing import normalize
import utils.tools as tl

base_path = '/Users/zhongyuanke/data/'
file_bulk = base_path + 'pbmc/bulk/GSE107011_Processed_data_TPM_symbol.csv'
seurat_celltype_path = base_path + 'multimodal/pbmc_10k/celltype_filt.csv'
file_vcca_label = base_path + 'vcca/atac/pbmc_10k/pvcca_01.h5ad'
file_orig_rna = base_path + 'multimodal/pbmc_10k/rna_filt.csv'
adata_bulk = sc.read_csv(file_bulk)
adata_vcca_label = sc.read_h5ad(file_vcca_label)
vcca_label = adata_vcca_label.obs['celltype']
seurat_label = tl.get_label_by_txt(seurat_celltype_path)
adata_orig = sc.read_csv(file_orig_rna)
print(adata_bulk, adata_orig)

# vcca_label = seurat_label

adata = adata_orig.concatenate(adata_bulk)
print(adata)
print(adata_bulk.obs)

x = adata.X[0:adata_orig.shape[0], ]
bulk_x = adata.X[adata_orig.shape[0]:adata.shape[0], ]
id_cMono = [i for i, x in enumerate(vcca_label) if x == 'Neutrophil']
cMono = x[id_cMono, ]
cMono = np.mean(cMono, axis=0)
print(cMono)

treg_bulk = bulk_x[26, ]

print(treg_bulk)

pccs = np.corrcoef(cMono, treg_bulk, )
print(pccs)

