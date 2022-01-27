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

resolution=2

# sc.pp.log1p(adata)
# sc.tl.leiden(adata, resolution=resolution)

rna.obsm['X_vcca']=adata.obsm['X_vcca']
sc.pp.log1p(rna)
sc.pp.neighbors(rna,use_rep='X_vcca', n_neighbors=13)
sc.tl.leiden(rna, resolution=resolution)

# dataframe = pd.DataFrame(adata.obs['leiden'])
# dataframe.to_csv("/Users/zhongyuanke/data/vcca/atac/pbmc_10k/vcca_leiden.csv",index=False,sep=',')
# print('finishi write leiden')
# sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
# sc.pl.rank_genes_groups(adata, n_genes=5, sharey=False,)
# --------------------------activity--------------------------------
# adata_atac = sc.read_h5ad(file_atac)
adata_atac = sc.read_csv(file_activity)
sc.pp.log1p(adata_atac)
print(adata_atac)
print(adata)
adata_atac.obsm['X_vcca'] = adata.obsm['X_vcca']
sc.pp.neighbors(adata_atac, use_rep='X_vcca', n_neighbors=13)
sc.tl.leiden(adata_atac, resolution=resolution)# resolution越大，簇就越多
#
# adata_atac.write_h5ad(base_path + 'vcca/atac/pbmc_10k/vcca_atac_temp.h5ad')
# sc.tl.rank_genes_groups(adata_atac, 'leiden', method='t-test')
# sc.pl.rank_genes_groups(adata_atac, n_genes=5, sharey=False,)
# sc.tl.umap(adata_atac)
# sc.pl.umap(adata_atac, color=['leiden'], s=7, legend_loc='on data')

# kmeans = KMeans(n_clusters=27, random_state=123, n_jobs=4).fit(adata.obsm['X_vcca'])
# kmeans_label = kmeans.labels_
# print(Counter(kmeans_label))
# kmeans_uns = pd.Series(kmeans_label, dtype='category')
# print(kmeans_uns)
# adata.obs['kmeans'] = kmeans_uns.values
# adata.uns['kmeans'] = kmeans_uns

print(adata)
# adata_test = sc.read_h5ad(file_vcca2)
# adata_protein = sc.read_csv(file_protein)

cluster2annotation_01_l2_14 = {
    '0':    'CD14 Mono',
    '1':	'CD14 Mono',
    '2':	'CD14 Mono',
    '3':    'CD8 TEM',
    '4':	'CD16 Mono',
    '5':	'CD8 Naive',
    '6':	'CD4 Naive',
    '7':	'CD14 Mono',
    '8':	'Memory B',
    '9':	'CD14 Mono',
    '10':	'CD8 TEM', #  *
    '11':	'CD8 Naive',
    '12':	'NK',
    '13':	'Treg',
    '14':	'CD14 Mono',
    '15':	'CD4 TCM',
    '16':	'CD4 TCM',
    '17':	'Naive B',
    '18':	'CD4 TEM',
    '19':	'CD4 Naive',
    '20':	'CD8 Naive',
    '21':	'CD8 Naive',
    '22':	'CD4 TCM',
    '23':	'CD4 Naive',
    '24':	'CD14 Mono',
    '25':	'CD4 TSCM', #*
    '26':	'CD14 Mono',
    '27':	'CD4 Naive',
    '28':	'CD8 Naive',
    '29':	'CD4 Naive',
    '30':	'cDC',
    '31':	'MAIT',
    '32':	'CD4 TCM',
    '33':	'Naive B',
    '34':	'pDC',
    '35':	'gdT',
    '36':	'NK',
    '37':	'CD14 Mono',
    '38':	'Naive B',
    '39':	'Basophil',
    '40':	'pre B',
    '41':	'HSPC',
    '42':   'Plasma'
}

# sc.tl.umap(adata)
# sc.pl.umap(adata, color=['CD11C','MHC','CD317','CD8A','CD4','CD11B'], s=2, ncols=4)
# sc.pl.umap(adata, color='leiden',legend_loc='on data')
# adata.write_h5ad(file_vcca)
# sc.pl.umap(adata, color=['seurat_celltype'], legend_loc='on data', frameon=False, s=5)

# sc.pl.umap(adata, color=['S100A8', 'CST3', 'CD8A', 'CD79A', 'CD74', 'MS4A1', 'LTB', 'FCER1G',
#                         'CD3E', 'CCL5', 'TYROBP', 'TCF7', 'LYZ', 'IL7R', 'HLA-DRA', 'LDHB'], s=2, ncols=4)
# # sc.pp.neighbors(rna)
# sc.tl.umap(rna)
# sc.pl.umap(rna, color=['TRPC3'], s=6, ncols=4)
# sc.pl.umap(adata, color=['TMEM74'], s=8, ncols=4)
# sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
# sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False)

# adata.obs['wnn_label'] = label
# adata.obs['vcca_celltype'] = adata.obs['kmeans'].map(cluster2annotation_kmeans_pvcca01).astype('category')


adata.obs['vcca_celltype_l2'] = adata.obs['leiden'].map(cluster2annotation_01_l2_14).astype('category')

celltype=pd.DataFrame(data=adata.obs['vcca_celltype_l2'].values)
celltype.to_csv('/Users/zhongyuanke/data/vcca/atac/pbmc_10k/vcca_celltype.csv', index=False,sep=',',encoding='gbk')

print('finishi write vcca_celltype')


sc.pl.umap(adata, color='vcca_celltype_l2', legend_loc='on data', frameon=False, s=8)
# adata.obs['vcca_celltype_l2'] = adata.obs['leiden'].map(cluster2annotation_01_l2).astype('category')
# rna.obs['vcca_celltype_l2'] = adata.obs['vcca_celltype_l2']
# rna.obs['vcca_leiden'] = adata.obs['leiden']
# rna.write_h5ad(file_rna)
adata.write_h5ad(file_vcca)
# adata_atac.obs['vcca_celltype_l2'] = adata_atac.obs['leiden'].map(cluster2annotation_01_l2).astype('category')
# adata_atac.obs['vcca_leiden'] = adata.obs['leiden']
# # sc.pl.umap(adata, color='vcca_celltype', legend_loc='on data', frameon=False, s=8)
# adata_atac.write_h5ad(file_atac)
# print(adata)
# adata_test.obs['celltype'] = adata.obs['celltype']
# sc.pl.umap(adata_test, color='celltype', legend_loc='on data', frameon=False, s=10)
#


# -------------------- write celltype ------------------------
# print(set(adata.obs['vcca_celltype_l2']))
# #字典中的key值即为csv中列名
# dataframe = pd.DataFrame(adata_atac.obs['vcca_celltype_l2'])
#
# #将DataFrame存储为csv,index表示是否显示行名，default=True
# dataframe.to_csv("/Users/zhongyuanke/data/vcca/atac/pbmc_10k/vcca_celltype.csv",index=False,sep=',')
#

marker_map_rna = {
    'neutrophil': ['LYST','CSF3R','MTRNR2L12'],
    'B cell':['CD19'],
    'B Memory': ['BANK1', 'CD79A', 'PLEKHG1','CD27'],
    'B Naive': ['BACH2', 'SYN3','CD38'],
    'CD27+CD38h B cells':['CD27','CD38','CD24','CD19'],
    'Plasma': ['CD44','CD38'],
    'Basophil': ['RAB31'],
    'CD14 Mono': ['CSF3R','VCAN','CD14','CCR2','CD163'],
    'CD16 Mono': ['FCGR3A','FCER1G', 'LYPD2','CX3CR1'],
    'Intermidiate Mono': ['CX3CR1','CCR5','FCGR3A','CSF1R','SELL','CCR2'],
    'temp Mono':['SKAP1','BCL11B','CD247'],
    'TCM':['CCR7','IL2RA','SELL','FOXP1','TCF7'],
    'TEM':['TBX21','GZMA'],
    'TEM1/TEM2':['CD27','CD28'],
    'naive T': ['CCR7'],
    'CD4 Naive': ['CD4','LEF1', 'CAMK4','FHIT', 'PTPRC','CCR7','GZMA'],
    'CD4 TEM': ['CD4','INPP4B','CD58'],
    'CD4 TEMRA':['CD4','PTPRC'],
    'Th1': ['CCR1','CXCR3'],
    'Th2': ['CCR4','CXCR4'],
    'CD4 TCM': ['CCR7','IL2RA','CD58','FOXP1','CD127','CD62L'],
    'CD8 Naive': ['CD8A','CD8B','CCR7','GZMA'],
    'CD8 TEM': ['CD8A','CCL5','CCR5','GZMA'],
    'CD8 TCM': ['CD8A','CCR7','CCR5','CD127','CD62L'],
    'MAIT': ['ME1', 'SLC4A10'],
    # 'CD27+CD38h B cells':['CD27','CD38','CD24','CD19'],
    'NK2': ['BCL11B'],
    'Treg2': ['CD27'],
    # 'Neutrophils': ['PLXDC2', 'CSF3R', 'RBM47','CXCL8'],
    # 'eosinophil':['SYNE2'],
    'NK T-cell': ['CACNA1C','CCR7','NCAM1','CD38'],
    # 'Myeloid DC': ['TPM2'],
    'Treg': ['STAM','CCR4','IL32','IL2RA'],
    'gdT': ['CST7', 'NKG7','GZMA'],
    'pDC': ['ZFAT', 'CD4'],
    'cDC': ['CD1C',  'CD33',],
    # 'iDC':[],
    # 'MDSC': ['CD15','CD33','CD66B'],

}
marker_map_atac = {
    'basophil': ['RAB31'],
    'pDC': ['PLD4','TXNDC5','SOX4','ZFAT','TCF4'],
    # 'neutrophil': ['MCTP2','PLXDC2','CSF3R','SULF2','MCTP2','MTRNR2L12'],
    'Non-classical Monocyte':['FCER1G', 'AOAH', 'CABP4', 'LYPD2'],
    'classical Monocyte':  ['FCN1', 'LYZ','S1PR3', 'RXFP2','VCAN','DMXL2'],

    'gdT': ['CST7','NKG7','JAKMIP1'],
    'Treg':  ['IL7R', 'IL32', 'STAM'],
    'CD8 T-cell': 'CD8B',
    'naive CD4 T-cell': ['BCL11B', 'LEF1'],
    'Memory CD4 T-cell': 'INPP4B',
    'naive B': ['BACH2','SYN3','CD74'],
    'NK T-cell': ['CACNA1C'],
    'MAIT': ['ME1','SLC4A10','ZBTB16'],
    'memory B-cell': ['PPP1R37','BANK1','CD79A','PLE    `KHG1'],
    'memory CD8 T-cell':['CCL5','CD8A'],
    'Myeloid DC':['TPM2'],
}

marker_map_rna_save = {
    # 'neutrophil': ['LYST','CSF3R','MTRNR2L12'],
    'B cell': ['CD19'],
    'B Memory': ['BANK1', 'CD79A', 'PLEKHG1'],
    'B Naive': ['BACH2', 'SYN3'],
    'B_pre': ['CD44', 'STAT3'],
    'Basophil': ['RAB31'],
    'CD14 Mono': ['CSF3R', 'VCAN', 'CD14','NEAT1'],
    'CD16 Mono': ['FCGR3A', 'FCER1G', 'LYPD2', 'CX3CR1'],
    # 'Inte Mono': ['CX3CR1', 'SELL','CCR2'],
    # 'temp Mono':['SKAP1','BCL11B','CD247'],
    # 'TCM':['CCR7','IL2RA','SELL','FOXP1','TCF7'],
    # 'TEM':['TBX21','GZMA'],
    # 'TEM1/TEM2':['CD27','CD28'],
    'CD4 T': ['CD4'],
    # 'Naive T': ['CCR7'],
    'CD4 Naive': ['LEF1', 'FHIT'],
    'CD4 TCM': ['CCR7', 'FOXP1'],
    'CD4 TEM': ['INPP4B',],
    # 'CD4 TEMRA':['CD4','PTPRC'],
    # 'Th1': ['CCR1','CXCR3'],
    # 'Th2': ['CCR4','CXCR4'],

    'CD4 TSCM': ['LEF1','CCR7'],
    'CD8 T': ['CD8A',],
    'CD8 Naive': ['CD8B', 'CCR7'],
    'CD8 TEM': ['CCL5', 'GZMA'],
    # 'CD8 TCM':['CCR7','CCR5'],
    'HSPC': ['CD34', 'RPS3A',],
    'MAIT': ['ME1', 'SLC4A10'],
    # 'CD27+CD38h B cells':['CD27','CD38','CD24','CD19'],
    # 'NK2':['BCL11B'],
    # 'Treg2':['CD27'],
    # 'Neutrophils': ['PLXDC2', 'CSF3R', 'RBM47','CXCL8'],
    # 'eosinophil':['SYNE2'],
    'NK T-cell': ['NCAM1'],
    # 'Myeloid DC': ['TPM2'],
    'Plasma': ['CD38'],
    'Treg': ['CCR4','IL32','IL2RA'],
    'cDC': ['CD1C',  'CD33',],
    'gdT': ['CST7', 'NKG7','GZMA'],
    'pDC': ['ZFAT', 'CD4'],
    # 'pre-B': ['CD44', 'STAT3'],
    # 'iDC':[],
    # 'MDSC': ['CD15','CD33','CD66B'],
}
# sc.pl.matrixplot(adata, marker_map, 'leiden', dendrogram=True, cmap='Blues',)
#
# sc.tl.umap(rna)
# sc.pl.umap(rna,color='leiden',legend_loc='on data')
rna.obs['vcca_celltype_l2'] = adata.obs['vcca_celltype_l2']
rna.layers['scaled'] = sc.pp.scale(rna, copy=True).X
sc.pl.matrixplot(rna, marker_map_rna_save, 'vcca_celltype_l2',
                 dendrogram=False, layer='scaled', vmin=-2, vmax=2, cmap='RdBu_r')
# sc.pl.matrixplot(rna, marker_map_rna_save, 'leiden', dendrogram=True,
#                  layer='scaled', vmin=-2, vmax=2, cmap='RdBu_r')
# sc.pl.matrixplot(adata, marker_map_rna, 'seurat_celltype', dendrogram=True,
#                  layer='scaled', vmin=-2, vmax=2, cmap='RdBu_r')

# adata_atac.obs['seurat_celltype']=adata.obs['seurat_celltype']
# adata_atac.layers['scaled'] = sc.pp.scale(adata_atac, copy=True).X
# # sc.pl.matrixplot(adata_atac, marker_map_atac, 'vcca_celltype_l2', dendrogram=True,
# #                  layer='scaled', vmin=-2, vmax=2, cmap='RdBu_r')
# sc.pl.matrixplot(adata_atac, marker_map_atac, 'seurat_celltype', dendrogram=True,
#                  layer='scaled', vmin=-2, vmax=2, cmap='RdBu_r')

# adata_atac.layers['scaled'] = sc.pp.scale(adata_atac, copy=True).X
# sc.pl.matrixplot(adata_atac, marker_map, 'vcca_celltype', dendrogram=False,
#                  layer='scaled', vmin=-2, vmax=2, cmap='RdBu_r')
# sc.pl.umap(adata, color=['SORCS1'],
#            s=2, ncols=4) # iMonocyte

# sc.pl.umap(adata, color=['HNF1B'],
#            s=2, ncols=4) # balas
# sc.pl.umap(adata, color=['TPM2'],
#            s=2, ncols=4) # myeloid
# sc.pl.umap(adata, color=['THBS4'],
#            s=2, ncols=4) # Eosinophils

# sc.pl.dotplot(adata, marker_map,
#               groupby='wnn_label', dendrogram=True)

# ax = sc.pl.heatmap(adata, marker_map, groupby='celltype', cmap='viridis', dendrogram=True)
# sc.tl.rank_genes_groups(adata, 'celltype', method='t-test')
# sc.pl.rank_genes_groups(adata, n_genes=5, sharey=False)
# sc.pl.umap(adata, color=['MT-ATP6', 'CD74', 'LYZ', 'IL32', 'CD3D', 'CST7', 'LEF1', 'FCER1G', 'PLD4', 'FTL'],
#            s=2, ncols=4)

# sc.pl.umap(adata, color=['CYP1B1'],
#            s=2, ncols=4)
# adata_protein.obs['celltype'] = adata.obs['celltype']
# sc.pp.neighbors(adata_protein)
# sc.tl.umap(adata_protein)
# sc.pl.umap(adata_protein, color='celltype', legend_loc='on data', frameon=False, s=5)
#
# sc.pp.neighbors(adata)
# sc.tl.umap(adata)
# sc.pl.umap(adata, color='celltype', legend_loc='on data', frameon=False, s=5)
# sc.pl.umap(adata, color=['S100A8', 'CST3', 'CD8A', 'CD79A', 'CD74', 'MS4A1', 'LTB', 'FCER1G',
#                        'CD3E', 'CCL5', 'TYROBP', 'TCF7', 'LYZ', 'IL7R', 'HLA-DRA', 'LDHB'], s=2, ncols=4)

