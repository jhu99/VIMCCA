import scanpy as sc
import numpy as np
import anndata
from sklearn.preprocessing import normalize
import utils.tools as tl
from collections import Counter

base_path = '/Users/zhongyuanke/data/'
file_bulk = base_path + 'pbmc/bulk/GSE107011_Processed_data_TPM_symbol.csv'
seurat_celltype_path = '/Users/zhongyuanke/data/multimodal/atac_pbmc_10k/seurat_celltype_filt.csv'
file_vcca_label = base_path + 'vcca/atac/pbmc_10k/vcca_temp.h5ad'
file_orig_rna = base_path + 'multimodal/atac_pbmc_10k/rna_filt.csv'
adata_bulk = sc.read_csv(file_bulk)
adata_vcca_label = sc.read_h5ad(file_vcca_label)
vcca_label = adata_vcca_label.obs['vcca_celltype_l2']
seurat_label = tl.get_label_by_txt(seurat_celltype_path)
adata_orig = sc.read_csv(file_orig_rna)
print(adata_bulk, adata_orig)

print(set(seurat_label))
print(set(vcca_label))
# vcca_label = seurat_label

adata = adata_orig.concatenate(adata_bulk)
print(adata)
print(adata_bulk.obs)

x = adata.X[0:adata_orig.shape[0], ]
bulk_x = adata.X[adata_orig.shape[0]:adata.shape[0], ]


# cellname = 'Memory B'
# bulk_v = bulk_x[47, ]

# cellname = 'Naive B'
# bulk_v = bulk_x[44, ]
#
# cellname = 'CD14 Mono'
# bulk_v = bulk_x[49, ]
#
# cellname = 'CD16 Mono'
# bulk_v = bulk_x[51, ]
#
# cellname = 'CD4 Naive'
# bulk_v = bulk_x[41, ]
#
# cellname = 'CD4 TEM'
# bulk_v = bulk_x[42, ]
#
# cellname = 'CD8 Naive'
# bulk_v = bulk_x[28, ]
#
# cellname = 'CD8 TEM_2'
# bulk_v = bulk_x[31, ]
#
# cellname = 'gdT'
# bulk_v = bulk_x[33, ]
#
# cellname = 'MAIT'
# bulk_v = bulk_x[32, ]
#
# cellname = 'NK'
# bulk_v = bulk_x[52, ]
#
# cellname = 'pDC'
# bulk_v = bulk_x[53, ]
#
# cellname = 'Plasma'
# bulk_v = bulk_x[48, ]
#
cellname = 'Treg'
bulk_v = bulk_x[36, ]


vcca_ids = [i for i, x in enumerate(vcca_label) if x == cellname]
vcca_v = x[vcca_ids, ]
vcca_v = np.mean(vcca_v, axis=0)
pccs = np.corrcoef(vcca_v, bulk_v)
print('VCCA '+cellname+' :', pccs)

wnn_ids = [i for i, x in enumerate(seurat_label) if x == cellname]
wnn_v = x[wnn_ids, ]
wnn_v = np.mean(wnn_v, axis=0)
pccs = np.corrcoef(wnn_v, bulk_v)
print('WNN '+cellname+' :', pccs)

