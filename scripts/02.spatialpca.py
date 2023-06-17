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
from importlib import resources
from spider import SPIDER
op = SPIDER()

arr = [
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
]

input1 =  int(sys.argv[1])-1
ds = arr[input1][0]
sample_name = arr[input1][1]
out_f = f'/home/lishiying/data6/SPIDER-paper/SPIDER-paper/datasets/{ds}/{sample_name}/'
idata = anndata.read_h5ad(f'{out_f}/idata.h5ad')
adata = anndata.read_h5ad(f'{out_f}/adata.h5ad')

if ds=='DFPLC':
    k = idata.uns['cell_meta'][adata.uns['cluster_key']].nunique()
    i = 7
    j = 50
    interface_id = np.where(idata.obs_names == idata.obs[idata.obs['label']=='WM_WM'].index[0])[0][0]
elif ds=='breast_cancer':
    k = sum(idata.uns['cell_meta'][idata.uns['cell_meta'][adata.uns['cluster_key']] != 'undetermined'][adata.uns['cluster_key']].value_counts() > 20)
    i = k+2
    j = 20
    interface_id = np.where(idata.obs_names == idata.obs[idata.obs['label']=='invasive cancer_invasive cancer'].index[0])[0][0]
    
    
count_f = f'{out_f}/idata_full_count.csv'
meta_f = f'{out_f}/idata_full_meta.csv'
idata.to_df().to_csv(count_f)
idata.obs[['row', 'col']].to_csv(meta_f)
with resources.path("spider.R_script", "run_spatialPCA.R") as pw_fn:
    os.system(str(f'/bin/bash -c "{R_path} -f /{pw_fn} {meta_f} {count_f} {k} {j} {out_f} {interface_id} {i}"'))