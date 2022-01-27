import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
from pyecharts.charts import Sankey
from pyecharts import options as opts
import scanpy as sc
from collections import Counter


base_path = '/Users/zhongyuanke/data/'
file_vcca = base_path + 'vcca/atac/pbmc_10k/vcca_temp.h5ad'

adata = sc.read_h5ad(file_vcca)

df_seurat = adata.obs['seurat_celltype'].values
df_vcca = adata.obs['vcca_celltype_l2'].values
print(df_seurat)
seurat_list = []
vcca_list = []

for n in df_seurat:
    seurat_list.append(n.replace(' ', '_')+'(Seurat4)')
for n in df_vcca:
    vcca_list.append(n.replace(' ', '_')+'(VIMCCA)')

nodes=[]
for name in set(seurat_list):
    nodes.append({'name':name})

for name in set(vcca_list):
    nodes.append({'name': name})
print(nodes)


def label2cluster(df_vcca, df_seurat):
    my_dict = {i: [] for i in df_vcca}
    for i, j in enumerate(df_seurat):
        # print(i,j)
        my_dict[df_vcca[i]].append(j)
    return my_dict


my_dict = label2cluster(vcca_list,seurat_list)
print(my_dict)
links=[]

for vcca in my_dict:
    print(vcca)
    seurat_name = my_dict[vcca]

    count = Counter(seurat_name)
    for c in count:
        temp_dict = {}
        temp_dict['source'] = vcca
        temp_dict['target']=c
        temp_dict['value']=count[c]
        links.append(temp_dict)
print(len(links))


sanky = Sankey(init_opts=opts.InitOpts(width="1400px", height="500px")) #设置图表的宽度和高度
sanky.add(
    "sankey",
    nodes,#读取节点
    links,#读取路径
    linestyle_opt=opts.LineStyleOpts(opacity=0.2, curve=0.5, color="source"),#设置线条样式
    label_opts=opts.LabelOpts(position="right"),#设置标签配置项
    node_align=str( "justify"),#设置节点对齐方式：right，left,justify(节点双端对齐)
)
# sanky.set_global_opts(title_opts=opts.TitleOpts(title="Sankey Annotation"))#表名
sanky.render('sankey.html')
# .render_notebook()

