#!/usr/bin/env /home/lishiying/.conda/envs/spider1/bin/python
import sys
import os
import numpy as np
import pandas as pd
import anndata
import scanpy as sc
R_path = 'source /etc/profile;module load GCC/11.2.0 OpenMPI/4.1.1 R/4.2.0 Anaconda3/2022.05 R-bundle-Bioconductor/3.15-R-4.2.0;R'
from importlib.metadata import version
version('spider-st')
from spider import SPIDER
op = SPIDER()

arr = [
 ['PDAC', 'PDAC_A'],
 ['DFPLC', '151673'],
 ['DFPLC', '151507'],
 ['DFPLC', '151508'],
 ['DFPLC', '151509'],
 ['DFPLC', '151510'],
 ['DFPLC', '151669'],
 ['DFPLC', '151670'],
 ['DFPLC', '151671'],
 ['DFPLC', '151672'],
 ['DFPLC', '151675'],
 ['DFPLC', '151676'],
 ['breast_cancer', 'A1'],
 ['breast_cancer', 'B1'],
 ['breast_cancer', 'C1'],
 ['breast_cancer', 'D1'],
 ['breast_cancer', 'E1'],
 ['breast_cancer', 'F1'],
 ['breast_cancer', 'G2'],
 ['breast_cancer', 'H1'],
 ['mouse_embryo', 'embryo1_2'],
 ['mouse_embryo', 'embryo1_5'],
 ['mouse_embryo', 'embryo2_2'],
 ['mouse_embryo', 'embryo2_5'],
 ['mouse_embryo', 'embryo3_2'],
 ['mouse_embryo', 'embryo3_5'],
 ['mouse_lung', 'D2'],
 ['mouse_lung', 'E4'],
 ['mouse_lung', 'E5'],
 ['slide-seq-v2', 'mousebrain'],
 ]

input1 =  int(sys.argv[1])-1
ds = arr[input1][0]
sample_name = arr[input1][1]
out_f = f'../input_datasets/{ds}/{sample_name}/'
print(ds, sample_name)
adata = anndata.read_h5ad(f'{out_f}/adata.h5ad')
print(adata)
no_spatalk = False
if len(adata) > 10000:
    no_spatalk=True
idata = op.prep(adata, out_f, R_path, cluster_key=adata.uns['cluster_key'], is_human=adata.uns['is_human'], coord_type=adata.uns['coord_type'], no_spatalk=no_spatalk)
idata, meta_idata = op.find_svi(idata, out_f, R_path, alpha=0.3)
idata.write_h5ad(f'{out_f}/idata.h5ad')
try:
    meta_idata.write_h5ad(f'{out_f}/meta_idata.h5ad')
except:
    pass
os.remove(f'{out_f}/adata_count.csv')
os.remove(f'{out_f}/adata_meta.csv')