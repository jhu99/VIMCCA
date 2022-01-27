import scanpy as sc
import utils.tools as tl
import umap
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

base_path = '/Users/zhongyuanke/data/'

file_rna = base_path + 'multimodal/protein/vdj_v1_hs_pbmc2_5gex_protein_filtered_feature_bc_matrix.h5'
file_protein = base_path + 'multimodal/protein/protein.csv'
file_vcca = base_path + 'vcca/pbmc_protein/vcca_temp_04.h5ad'
# celltype_path = base_path + 'multimodal/atac_pbmc_10k/celltype_filt.csv'

# label = tl.get_label_by_txt(celltype_path)
# label = tl.text_label_to_number(label)

# file_rna = base_path + 'multimodal/protein/rna.csv'
# file_protein = base_path + 'multimodal/protein/protein.csv'
# file_vcca = base_path + 'vcca/pbmc_protein/pvcca_02.h5ad'
# file_vcca2 = base_path + 'vcca/pbmc_protein/pvcca_03.h5ad'


adata = sc.read_h5ad(file_vcca)
rna = sc.read_10x_h5(file_rna)
adata.var_names_make_unique()
# dataframe = pd.DataFrame(adata.obs['leiden'])
# dataframe.to_csv("/Users/zhongyuanke/data/vcca/pbmc_protein/vcca_leiden.csv",index=False,sep=',')
# print('finishi write leiden')
#
# sc.tl.louvain(adata, resolution=2)
# dataframe = pd.DataFrame(adata.obs['louvain'])
# dataframe.to_csv("/Users/zhongyuanke/data/vcca/pbmc_protein/louvain_leiden.csv",index=False,sep=',')
# print('finishi write louvain')

resolution=4
# sc.pp.log1p(adata)
sc.tl.leiden(adata, resolution=resolution)
# sc.pl.umap(adata,color='leiden',legend_loc='on data',s=2)
print(adata)
# adata_test = sc.read_h5ad(file_vcca2)
adata_protein = sc.read_csv(file_protein)


cluster2annotation_02_l1 = {
    '0': 'B',
    '1': 'Mono',
    '2': 'Mono',
    '3': 'CD8 T',
    '4': 'CD4 T',
    '5': 'B',
    '6': 'CD8 T',
    '7': 'Mono',
    '8': 'Mono',
    '9': 'NK',
    '10': 'CD4 T',
    '11': 'Mono',
    '12': 'CD4 T',
    '13': 'CD4 T',
    '14': 'Mono',
    '15': 'Mono',
    '16': 'CD4 T',
    '17': 'CD8 T',
    '18': 'Mono',
    '19': 'Mono',
    '20': 'CD4 T',
    '21': 'Mono',
    '22': 'CD8 T',
    '23': 'CD4 T',
    '24': 'CD4 T',
    '25': 'CD8 T',
    '26': 'B',
    '27': 'NK',
    '28': 'CD4 T',
    '29': 'CD4 T',
    '30': 'Mono',
    '31': 'Mono',
    '32': 'B',
    '33': 'Mono',
    '34': 'CD4 T',
    '35': 'CD4 T',
    '36': 'CD8 T',
    '37': 'Mono',
    '38': 'B',
    '39': 'Mono',
    '40': 'CD4 T',
    '41': 'Mono',
    '42': 'Mono',
    '43': 'CD8 T',
    '44': 'B',
    '45': 'CD4 T',
    '46': 'CD4 T',
    '47': 'CD4 T',
    '48': 'B',
    '49': 'Mono',
    '50': 'Mono',
    '51': 'CD4 T',
    '52': 'Mono',
    '53': 'Mono',
    '54': 'NK',
    '55': 'B',
    '56': 'Mono',
    '57': 'Mono',
    '58': 'B',
    '59': 'Plasma',
    '60': 'CD4 T',
    '61': 'B',
    '62': 'B',
    '63': 'NK',
}


cluster2annotation_02_l2 = {
    '0': 'B Naive',
    '1': 'CD14 Mono',
    '2': 'CD14 Mono',
    '3': 'CD8 Naive',
    '4': 'CD4 Naive',
    '5': 'B Naive',
    '6': 'MAIT',
    '7': 'CD16 Mono',
    '8': 'CD14 Mono',
    '9': 'NK',
    '10': 'CD4 TEM(RO+)',
    '11': 'CD14 Mono',
    '12': 'CD4 Naive',
    '13': 'CD4 Naive',
    '14': 'CD14 Mono',
    '15': 'CD14 Mono',
    '16': 'CD4 Naive',
    '17': 'CD8 Naive',
    '18': 'CD14 Mono',
    '19': 'CD14 Mono',
    '20': 'CD4 Naive',
    '21': 'CD14 Mono',
    '22': 'CD8 TEM',
    '23': 'CD4 Naive',
    '24': 'CD4 Naive',
    '25': 'gdT',
    '26': 'B Naive',
    '27': 'NK',
    '28': 'CD4 Naive',
    '29': 'CD4 Naive',
    '30': 'CD14 Mono',
    '31': 'CD14 Mono',
    '32': 'B Memory',
    '33': 'CD14 Mono',
    '34': 'CD4 Naive',
    '35': 'CD4 TEM(RO+)',
    '36': 'CD8 Naive',
    '37': 'CD14 Mono',
    '38': 'B Naive',
    '39': 'CD14 Mono',
    '40': 'Treg',
    '41': 'cDC',
    '42': 'CD14 Mono',
    '43': 'MAIT',
    '44': 'B Naive',
    '45': 'CD4 TEM(RA+)',
    '46': 'CD4 Naive',
    '47': 'CD4 Naive',
    '48': 'B Memory',
    '49': 'CD14 Mono',
    '50': 'CD14 Mono',
    '51': 'CD4 Naive',
    '52': 'Inte Mono',
    '53': 'CD14 Mono',
    '54': 'NK',
    '55': 'B Naive',
    '56': 'CD14 Mono(IgG)',
    '57': 'pDC',
    '58': 'B_CD3+',
    '59': 'Basophil',
    '60': 'CD4 Naive',
    '61': 'CD14 Mono(RA+)',
    '62': 'Plasma',
    '63': 'NK',
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
# sc.pl.rank_genes_groups(adata, n_genes=8, sharey=True,ncols=8)
# adata_protein.obs['leiden'] = adata.obs['leiden']
# sc.tl.rank_genes_groups(adata_protein, 'leiden', method='t-test')
# sc.pl.rank_genes_groups(adata_protein, n_genes=4, sharey=True, ncols=10)
# adata.obs['wnn_label'] = label
# adata.obs['vcca_celltype_l1'] = adata.obs['leiden'].map(cluster2annotation_02_l1).astype('category')
# sc.pl.umap(adata, color='vcca_celltype_l1', legend_loc='on data', frameon=False, s=7)
adata.obs['vcca_celltype_l2'] = adata.obs['leiden'].map(cluster2annotation_02_l2).astype('category')
sc.pl.umap(adata, color='vcca_celltype_l2', legend_loc='on data', frameon=False, s=7)
# adata.write_h5ad(file_vcca)
# print(adata)
# adata_test.obs['celltype'] = adata.obs['celltype']
# sc.pl.umap(adata_test, color='celltype', legend_loc='on data', frameon=False, s=10)

dataframe = pd.DataFrame(adata.obs['vcca_celltype_l2'])
dataframe.to_csv("/Users/zhongyuanke/data/vcca/pbmc_protein/vcca_celltype.csv",index=False,sep=',')
print('finishi write vcca_celltype')


marker_map_rna = {
    'eosinophil':['SYNE2'],
    'neutrophil': ['PROK2','PLXDC2', 'CSF3R', 'TREM1','SLC2A3'],
    'B cell':['CD19','CD79A','MS4A1','IGHM','FAM129C','BANK1', 'PLEKHG1'],
    'B Memory': ['CD80','TAS1R3','CXCR3','CD27','CD40',],
    'B cells subtype':['CD27','CD38','CD24','CD19','CD74'],
    'Regular B':['CD1D','CD38'],
    'B Naive': ['BACH2', 'SYN3','CD38'],
    'Plasma': ['CD44'],
    'Basophil': ['RAB31','PPBP','MT-ATP6','SLC2A3','PPP1R10'],
    'CD14 Mono': ['GAS7','VCAN','CD14','CCR2','CD163','S100A8','S100A9','LYZ'],
    'CD16 Mono': ['FCER1G', 'LYPD2','CX3CR1','FCGR3A'],
    'Inte Mono(CD14++CD16+)': ['UACA','CX3CR1','CCR5','FCGR3A','CSF1R','SELL','CCR2'],
    'naive T': ['CCR7','LEF1'],
    'TEM': ['GZMA',],
    'TCM':['CD58','CD69','IL2RA','CD27'],
    'TSCM':['SELL','FOXP1','LEF1'],
    'CD4 T':['CD4'],
    'CD4 Naive': ['LEF1', 'CAMK4','FHIT','ABLIM1'],
    'CD4 TEM': ['INPP4B','CCR5','GZMA','CD58'],
    'Th1': ['CCR1','CXCR3'],
    'Th2': ['CCR4','CXCR4'],
    'CD4 TCM': ['CCR7','CCR5','IL2RA','CD58','FOXP1'],
    'CD8 T':['CD8A','CD8B'],
    'CD8 Naive': ['CD8B','CCR7'],
    'CD8 TEM': ['CCL5','CCR5'],
    'CD8 TCM':['CCR7','CCR5'],
    'MAIT': ['ME1', 'SLC4A10','TC2N'],
    # 'Neutrophils': ['PLXDC2', 'CSF3R', 'RBM47','CXCL8'],
    # 'eosinophil':['SYNE2'],
    'NK': ['CACNA1C','CCR7','NCAM1','CD38','KLRD1'],
    # 'Myeloid DC': ['TPM2'],
    'Treg': ['STAM','CD28','SIRPG','IL2RA','FOXP3','CD27','IL32'],
    'gdT': ['CST7', 'NKG7','CCL5','KLRD1','S100B'],
    'pDC': ['ZFAT', 'TCF4','CD4','FAM129C','STMN1'],
    'cDC': ['CD1C',  'CD33', 'CD2','CD40'],
    'cDC1': ['CD8A'],
    'myeloid DC':['LYZ'],
    # 'cDC2': ['CD11B'],
    # 'iDC':[],
    # 'MDSC': ['CD15','CD33','CD66B'],
    'MDC1': [ 'CD33']
}

marker_map_rna_save = {
    # 'eosinophil':['SYNE2'],
    # 'neutrophil': ['PROK2','PLXDC2', 'CSF3R', 'TREM1','SLC2A3'],
    'B cell':['CD19','CD79A','MS4A1','IGHM','FAM129C','BANK1', 'PLEKHG1'],
    'B Memory': ['CD80','CD27','CD40',],
    # 'B cells subtype':['CD27','CD38','CD24','CD19','CD74'],
    # 'Regular B':['CD1D','CD38'],
    'B Naive': ['BACH2', 'SYN3','CD38'],
    'Basophil': ['RAB31','PPBP'],
    'CD14 Mono': ['GAS7','VCAN','CD163','S100A8','S100A9','LYZ'],
    'CD16 Mono': ['FCER1G', 'LYPD2','CX3CR1','FCGR3A'],
    # 'Inte Mono(CD14++CD16+)': ['UACA','CX3CR1','CCR5','FCGR3A','CSF1R','SELL','CCR2'],
    'naive T': ['CCR7','LEF1'],
    # 'TCM':['CD58','CD69','IL2RA','CD27'],
    # 'TSCM':['SELL','FOXP1','LEF1'],
    # 'CD4 T':['CD4'],
    'CD4 Naive': ['LEF1', 'FHIT', 'ABLIM1'],
    'CD4 TEM': ['INPP4B'],
    # 'Th1': ['CCR1','CXCR3'],
    # 'Th2': ['CCR4','CXCR4'],
    'CD4 TCM': ['CCR7','CCR5','IL2RA','CD58','FOXP1'],
    'CD8 T':['CD8A','CD8B'],
    'CD8 Naive': ['CD8B','CCR7'],
    'CD8 TEM': ['CCL5','CCR5'],
    'CD8 TCM':['CCR7','CCR5'],
    'MAIT': ['ME1', 'SLC4A10','TC2N'],
    # 'Neutrophils': ['PLXDC2', 'CSF3R', 'RBM47','CXCL8'],
    # 'eosinophil':['SYNE2'],
    'NK': ['CD38','KLRD1'],
    # 'Myeloid DC': ['TPM2'],
    'Plasma': ['CD38','CD27'],
    'Treg': ['STAM','CD28','SIRPG','IL2RA','FOXP3','IL32'],
    'cDC': ['CD1C',  'CD33', ],
    'gdT': ['CST7', 'NKG7','CCL5','KLRD1','S100B'],
    'pDC': ['ZFAT', 'TCF4','CD4','FAM129C','STMN1'],
    # 'cDC1': ['CD8A'],
    # 'myeloid DC':['LYZ'],
    # 'cDC2': ['CD11B'],
    # 'iDC':[],
    # 'MDSC': ['CD15','CD33','CD66B'],
    # 'MDC1': [ 'CD33']
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

marker_map_protein = {'Protein':['CD25_TotalSeqC','CD14_TotalSeqC','CD16_TotalSeqC','CD19_TotalSeqC',
                      'IgG1_control_TotalSeqC','IgG2a_control_TotalSeqC','IgG2b_control_TotalSeqC','TIGIT_TotalSeqC',
                      'CD127_TotalSeqC','CD3_TotalSeqC','CD45RA_TotalSeqC','CD45RO_TotalSeqC',
                      'CD8a_TotalSeqC', 'PD-1_TotalSeqC','CD25_TotalSeqC',
                      'CD56_TotalSeqC','CD15_TotalSeqC','CD4_TotalSeqC']
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
adata.layers['scaled'] = sc.pp.scale(adata, copy=True).X
# sc.pl.matrixplot(adata, marker_map_rna_save, 'leiden', dendrogram=False,
#                  layer='scaled', vmin=-2, vmax=2, cmap='RdBu_r')
# sc.pl.matrixplot(adata, marker_map_rna_save, 'vcca_celltype_l2', dendrogram=False,
#                  layer='scaled', vmin=-2, vmax=2, cmap='RdBu_r')

adata_protein.obsm['X_vcca'] = adata.obsm['X_vcca']
sc.pp.neighbors(adata_protein, use_rep='X_vcca', n_neighbors=13)
sc.tl.leiden(adata_protein, resolution=resolution)
# sc.tl.umap(adata_protein)
print(adata_protein.var_names)
adata_protein.obs['vcca_celltype_l2'] = adata.obs['vcca_celltype_l2']
# adata_protein.write_h5ad(base_path+'vcca/pbmc_protein/protein_leiden.h5ad')
adata_protein.layers['scaled'] = sc.pp.scale(adata_protein, copy=True).X
sc.pl.matrixplot(adata_protein, marker_map_protein, 'vcca_celltype_l2', dendrogram=False,
                 layer='scaled', vmin=-2, vmax=2, cmap='RdBu_r')
# sc.pl.matrixplot(adata_protein, marker_map_protein, 'leiden', dendrogram=True, cmap='RdBu_r',layer='scaled',
#                  vmin=-2, vmax=2)
