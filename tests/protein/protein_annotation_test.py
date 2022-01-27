import scanpy as sc
import utils.tools as tl
import umap
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

base_path = '/Users/zhongyuanke/data/'

file_rna = base_path + 'multimodal/protein/vdj_v1_hs_pbmc2_5gex_protein_filtered_feature_bc_matrix.h5'
file_protein = base_path + 'multimodal/protein/protein.csv'
file_vcca = base_path + 'vcca/pbmc_protein/vcca_temp_01.h5ad'
# celltype_path = base_path + 'multimodal/atac_pbmc_10k/celltype_filt.csv'

# label = tl.get_label_by_txt(celltype_path)
# label = tl.text_label_to_number(label)

# file_rna = base_path + 'multimodal/protein/rna.csv'
# file_protein = base_path + 'multimodal/protein/protein.csv'
# file_vcca = base_path + 'vcca/pbmc_protein/pvcca_02.h5ad'
# file_vcca2 = base_path + 'vcca/pbmc_protein/pvcca_03.h5ad'


adata = sc.read_h5ad(file_vcca)
rna = sc.read_10x_h5(file_rna)

# dataframe = pd.DataFrame(adata.obs['leiden'])
# dataframe.to_csv("/Users/zhongyuanke/data/vcca/pbmc_protein/vcca_leiden.csv",index=False,sep=',')
# print('finishi write leiden')
#
# sc.tl.louvain(adata, resolution=2)
# dataframe = pd.DataFrame(adata.obs['louvain'])
# dataframe.to_csv("/Users/zhongyuanke/data/vcca/pbmc_protein/louvain_leiden.csv",index=False,sep=',')
# print('finishi write louvain')

resolution=2
# sc.pp.log1p(adata)
# sc.tl.leiden(adata, resolution=resolution)
# sc.pl.umap(adata,color='leiden',legend_loc='on data')
print(adata)
# adata_test = sc.read_h5ad(file_vcca2)
adata_protein = sc.read_csv(file_protein)

cluster2annotation_02_l2 = {
    '0': 'B Naive',
    '1': 'CD4 TEM',
    '2': 'CD14 Mono',
    '3': 'CD14 Mono',
    '4': 'CD4 TEM',
    '5': 'MAIT',
    '6': 'CD4 Naive',
    '7': 'NK',
    '8': 'CD4 Naive',
    '9': 'CD14 Mono',
    '10': 'CD14 Mono',
    '11': 'CD8 TEM',
    '12': 'Inte Mono',
    '13': 'CD14 Mono',
    '14': 'B Naive',
    '15': 'CD14 Mono',
    '16': 'CD4 Naive',
    '17': 'CD14 Mono(igG1)',
    '18': 'B Memory',
    '19': 'CD8 Naive',
    '20': 'CD4 Naive',
    '21': 'CD16 Mono',
    '22': 'Treg',
    '23': 'CD4 Naive',
    '24': 'CD4 Naive',
    '25': 'CD14 Mono',
    '26': 'B Naive',
    '27': 'cDC',
    '28': 'pDC',
    '29': 'Basophil',
    '30': 'Plasma',

}

# sc.tl.leiden(adata)
# adata.write_h5ad(file_vcca)
# sc.pl.umap(adata, color='leiden', legend_loc='on data', frameon=False, s=5)
# sc.pl.umap(adata, color=['S100A8', 'CST3', 'CD8A', 'CD79A', 'CD74', 'MS4A1', 'LTB', 'FCER1G',
#                         'CD3E', 'CCL5', 'TYROBP', 'TCF7', 'LYZ', 'IL7R', 'HLA-DRA', 'LDHB'], s=2, ncols=4)
# sc.pp.neighbors(rna)
# sc.tl.umap(rna)
# sc.pl.umap(rna, color=['TRPC3'], s=6, ncols=4)
# sc.pl.umap(adata, color=['TMEM74'], s=8, ncols=4)
# sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
# sc.pl.rank_genes_groups(adata, n_genes=8, sharey=True)
adata_protein.obs['leiden'] = adata.obs['leiden']
# sc.tl.rank_genes_groups(adata_protein, 'leiden', method='t-test')
# sc.pl.rank_genes_groups(adata_protein, n_genes=4, sharey=True)
# adata.obs['wnn_label'] = label
adata.obs['vcca_celltype'] = adata.obs['leiden'].map(cluster2annotation_02_l2).astype('category')
# sc.pl.umap(adata, color='vcca_celltype', legend_loc='on data', frameon=False, s=5)
# adata.write_h5ad(file_vcca)
# print(adata)
# adata_test.obs['celltype'] = adata.obs['celltype']
# sc.pl.umap(adata_test, color='celltype', legend_loc='on data', frameon=False, s=10)
#
# marker_map = {
#     'pDC': ['PLD4'],
#     'Non-classical Monocyte': ['FCER1G', 'AOAH', 'CABP4', 'LYPD2'],
#     'classical Monocyte': ['FCN1', 'S1PR3', 'RXFP2'],
#     'B cell': ['CD79A'],
#     'gdT': 'CST7',
#     'Treg':  ['IL7R', 'IL32', ],
#     'CD8 T-cell': 'CD8B',
#     'naive CD4 T-cell': 'FHIT',
#     'Memory CD4 T-cell': 'INPP4B',
#     'naive B': ['BACH2','SYN3'],
#     'NK T-cell': ['CD8A', 'CD8B'],
#     'MAIT': ['ME1','SLC4A10'],
#     'memory B': 'PPP1R37',
# }

marker_map_rna = {
    'basophil': ['RAB31'],
    'classical Monocyte': ['FCN1', 'GAS7', 'VCAN', 'DMXL2'],
    'MAIT': ['ME1', 'SLC4A10','FKBP11'],
    'memory B-cell': ['BANK1', 'CD79A', 'PLEKHG1'],
    'Memory CD4 T-cell': 'INPP4B',
    'memory CD8 T-cell': ['CCL5', 'PRKCH'],
    'naive B': ['BACH2', ],
    'naive CD4 T-cell': ['LEF1', 'CAMK4','FHIT'],
    'naive CD8 T-cell': ['CD8B'],
    'neutrophil': ['PLXDC2', 'CSF3R', 'TREM1','IQCN'],
    'Non-classical Monocyte': ['FCER1G', 'LYPD2', 'CSF1R'],
    'Treg': ['STAM', 'CD28', 'SIRPG'],
    # 'eosinophil':['SYNE2'],
    'NK T-cell': ['CACNA1C', 'KLRD1'],
    # 'Myeloid DC': ['TPM2'],
    'gdT': ['CST7', 'NKG7', ],
    'pDC': ['ZFAT', 'TCF4',],
}
marker_map_rna2 = {
    'B Memory': ['BANK1', 'CD79A', 'PLEKHG1','CD27'],
    'B Naive': ['BACH2', 'SYN3','CD19'],
    'Basophil': ['RAB31'],
    'C Mono': ['GAS7','VCAN'],
    'CD4 Naive': ['BCL11B', 'LEF1', 'CAMK4','FHIT'],
    'CD4 TEM': 'INPP4B',
    'CD8 Naive': ['CD8B'],
    'CD8 TEM': ['CCL5','CCR7',],
    'MAIT': ['ME1', 'SLC4A10'],
    'NC Mono': ['FCER1G', 'LYPD2'],
    # 'Neutrophils': ['PLXDC2', 'CSF3R', 'RBM47','CXCL8'],
    # 'eosinophil':['SYNE2'],
    'NK T-cell': ['CACNA1C','CCR7'],
    # 'Myeloid DC': ['TPM2'],
    'Treg': ['STAM'],
    'gdT': ['CST7', 'NKG7'],
    'pDC': ['ZFAT', 'TCF4','CD4'],
    'cDC': ['CD1C',  'CD33', 'CD2', ],
    'MDC1': ['CD33']
}
marker_map_rna3 = {
    'neutrophil': ['PROK2','PLXDC2', 'CSF3R', 'TREM1',],
    'B cell':['CD19'],
    'B Memory': ['BANK1', 'CD79A', 'PLEKHG1','CD27','CD80'],
    'Regular B':['CD1D'],
    'B Naive': ['BACH2', 'SYN3','CD38'],
    'Plasma': ['CD44',],
    'Basophil': ['RAB31','PPBP'],
    'CD14 Mono': ['GAS7','VCAN','CD14','CCR2','CD163'],
    'CD16 Mono': ['FCER1G', 'LYPD2','CX3CR1','FCGR3A'],
    'Inte Mono(CD14++CD16+)': ['UACA','CX3CR1','CCR5','FCGR3A','CSF1R','SELL','CCR2'],
    'naive T': ['CCR7'],
    'CD4 Naive': ['LEF1', 'CAMK4','FHIT'],
    'CD4 TEM': ['INPP4B','CCR5','GZMA','CD58'],
    'Th1': ['CCR1','CXCR3'],
    'Th2': ['CCR4','CXCR4'],
    'CD4 TCM': ['CCR7','CCR5','IL2RA','CD58','FOXP1'],
    'CD8 Naive': ['CD8B','CCR7'],
    'CD8 TEM': ['CCL5','CCR5'],
    'CD8 TCM':['CCR7','CCR5'],
    'MAIT': ['ME1', 'SLC4A10'],
    'CD27+CD38h B cells':['CD27','CD38','CD24','CD19'],
    # 'Neutrophils': ['PLXDC2', 'CSF3R', 'RBM47','CXCL8'],
    # 'eosinophil':['SYNE2'],
    'NK': ['CACNA1C','CCR7','NCAM1','CD38','CCL5','KLRD1'],
    # 'Myeloid DC': ['TPM2'],
    'Treg': ['STAM','CD28','SIRPG','IL2RA','FOXP3'],
    'gdT': ['CST7', 'NKG7','CCL5','TSPOAP1','KLRD1'],
    'pDC': ['ZFAT', 'TCF4','CD4'],
    'cDC': ['CD1C',  'CD33', 'CD2','CD40'],
    # 'iDC':[],
    # 'MDSC': ['CD15','CD33','CD66B'],
    'MDC1': [ 'CD33']
}
marker_map_protein = {
    'basophil': ['CD25_TotalSeqC'],
    'CD14 Mono': ['CD14_TotalSeqC'],
    'CD16 Mono': ['CD16_TotalSeqC'],
    'B memory': ['CD19_TotalSeqC','IgG1_control_TotalSeqC','IgG2a_control_TotalSeqC','IgG2b_control_TotalSeqC'],
    'B Naive': ['CD19_TotalSeqC', 'IgG1_control_TotalSeqC'],
    'T cell':['IgG2b_control_TotalSeqC','TIGIT_TotalSeqC'],
    'MAIT': ['CD127_TotalSeqC','CD3_TotalSeqC'],
    'CD4 Naive': ['CD3_TotalSeqC','CD45RA_TotalSeqC'],
    'CD4 TEM': ['CD127_TotalSeqC','CD3_TotalSeqC','CD45RO_TotalSeqC'],
    'CD8 TEM': ['CD8a_TotalSeqC','CD3_TotalSeqC', 'PD-1_TotalSeqC','CD45RO_TotalSeqC'],
    'CD8 Naive': ['CD3_TotalSeqC','CD45RA_TotalSeqC','PD-1_TotalSeqC'],
    # 'neutrophil': ['PLXDC2', 'CSF3R', 'TREM1','IQCN'],
    # 'Non-classical Monocyte': ['FCER1G', 'LYPD2', 'CSF1R'],
    'Treg': ['CD3_TotalSeqC','TIGIT_TotalSeqC','CD25_TotalSeqC'],
    # 'eosinophil':['SYNE2'],
    'NK': ['CD3_TotalSeqC','CD16_TotalSeqC','CD56_TotalSeqC','TIGIT_TotalSeqC'],
    'Myeloid DC': ['CD15_TotalSeqC'],
    'gdT': ['CD3_TotalSeqC'],
    'pDC': ['CD4_TotalSeqC'],
}
# sc.pp.filter_genes(rna, min_cells=1)
# sc.pp.log1p(rna)
# # sc.pp.scale(rna)
# rna.var_names_make_unique()
# rna.obsm['X_vcca']=adata.obsm['X_vcca']
# sc.pp.neighbors(rna, use_rep='X_vcca')
# sc.tl.leiden(rna, resolution=resolution)

# sc.tl.umap(rna)
# rna.obs['vcca_celltype']=adata.obs['vcca_celltype']
# rna.write_h5ad(base_path+'vcca/pbmc_protein/rna_leiden.h5ad')
# print(rna)
adata.var_names_make_unique()
adata.layers['scaled'] = sc.pp.scale(adata, copy=True).X
sc.pl.matrixplot(adata, marker_map_rna3, 'leiden', dendrogram=True,
                 layer='scaled', vmin=-2, vmax=2, cmap='RdBu_r')

adata_protein.var_names_make_unique()
# print(adata_protein.var.values)
adata_protein.obsm['X_vcca'] = adata.obsm['X_vcca']
sc.pp.neighbors(adata_protein, use_rep='X_vcca')
sc.tl.leiden(adata_protein, resolution=resolution)
# sc.tl.umap(adata_protein)
print(adata_protein.var_names)
adata_protein.obs['vcca_celltype'] = adata.obs['vcca_celltype']
# adata_protein.write_h5ad(base_path+'vcca/pbmc_protein/protein_leiden.h5ad')
adata_protein.layers['scaled'] = sc.pp.scale(adata_protein, copy=True).X
# sc.pl.matrixplot(adata_protein, marker_map_protein, 'vcca_celltype', dendrogram=False,
#                  layer='scaled', vmin=-2, vmax=2, cmap='RdBu_r')
sc.pl.matrixplot(adata_protein, marker_map_protein, 'leiden', dendrogram=True, cmap='RdBu_r',layer='scaled',
                 vmin=-2, vmax=2)
