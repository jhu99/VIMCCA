import scanpy as sc
from scbean.model import vcca as vcca
from scbean.tools import utils as tl
import matplotlib
import pandas as pd
matplotlib.use('TkAgg')

base_path = '/Users/zhongyuanke/data/'
file_rna = base_path + 'multimodal/protein/rna.csv'
file_rna = base_path + 'multimodal/protein/vdj_v1_hs_pbmc2_5gex_protein_filtered_feature_bc_matrix.h5'
file_protein = base_path + 'multimodal/protein/protein.csv'

adata_x= sc.read_10x_h5(file_rna)
adata_y = sc.read_csv(file_protein)


sc.pp.filter_genes(adata_x, min_cells=1)
sc.pp.log1p(adata_x)
sc.pp.log1p(adata_y)
sc.pp.scale(adata_x)
sc.pp.scale(adata_y)
# sc.pp.highly_variable_genes(adata_x, n_top_genes=5000)
# adata_x = adata_x[:, adata_x.var.highly_variable]

from concurrent.futures import ThreadPoolExecutor

with ThreadPoolExecutor() as executor:
    monitor = MemoryMonitor()
    mem_thread = executor.submit(monitor.measure_usage)
    try:
        fn_thread = executor.submit(desc.train, adata, dims=[adata.shape[1], 32, 16], tol=0.005, n_neighbors=10,
                   batch_size=256, louvain_resolution=[0.8],
                   save_dir="result", do_tsne=False, learning_rate=300,
                   do_umap=False,
                   save_encoder_weights=False)
        result = fn_thread.result()

    finally:
        monitor.keep_measuring = False
        max_usage = mem_thread.result()

    print(f"Peak memory usage: {max_usage/1024/1024/1024} GB")

z = vcca.fit_integration(
    adata_x, adata_y,
    sparse_x=False,
    sparse_y=False,
    mode='VCCA',
    hidden_layers=[128, 64, 32, 8],
    epochs=50
)

