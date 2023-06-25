import matplotlib.pyplot as plt
from spider import SPIDER
op = SPIDER()
import anndata
import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np

import os
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
TF_ENABLE_ONEDNN_OPTS = 0

import importlib
importlib.metadata.version('spider-st')


ds = 'mouse_embryo'
sample_name = 'embryo1_2'
out_f = f'../input_datasets/{ds}/{sample_name}/'
R_path = 'source /etc/profile;module load GCC/11.2.0 OpenMPI/4.1.1 R/4.2.0 Anaconda3/2022.05 R-bundle-Bioconductor/3.15-R-4.2.0;R'
adata = anndata.read_h5ad(f'{out_f}/adata.h5ad')
no_spatalk = False
if len(adata) > 10000:
    no_spatalk=True
idata = op.prep(adata, out_f, R_path, cluster_key=adata.uns['cluster_key'], is_human=adata.uns['is_human'], coord_type=adata.uns['coord_type'], no_spatalk=no_spatalk)
idata, meta_idata = op.find_svi(idata, out_f, R_path, alpha=0.3)
idata.write_h5ad(f'{out_f}/idata.h5ad')
idata = anndata.read_h5ad(f'{out_f}/idata.h5ad')


svi_df, svi_df_strict = op.svi.combine_SVI(idata,threshold=0.01)
op.svi.eva_SVI(idata, svi_df_strict)
plt.savefig(f'../figures/{ds}_{sample_name}_metric(supp).png', dpi=300,bbox_inches='tight')
plt.close()

op.util.adata_moranI(adata, out_f=out_f)
plt.rcParams['font.size'] = 13
merged_df,lri_pw_list,gene_lr_list,gene_pw_list = op.vis.svg_svi_relation(adata, idata, title='embryo 1 (z2)', is_human=adata.uns['is_human'], top=50)
plt.savefig(f'../figures/{ds}_{sample_name}_relation(main4D).png', dpi=600,bbox_inches='tight')
plt.close()


import gseapy
membership=pd.get_dummies(merged_df.set_index('Term')['group']).groupby('Term').sum().astype(str).agg('-'.join, axis=1).reset_index()
for x in membership[0].unique():
    sub_df = merged_df[merged_df.Term.isin(membership[membership[0]==x].Term)]
    if len(sub_df.Term.unique())>20:
        ax = gseapy.dotplot(sub_df, title='',cmap='viridis_r', size=20, figsize=(3,8),top_term=30,show_ring=True, ofname=f'../figures/{ds}_{sample_name}_svgi_pw_{x}(supp).png')
    elif len(sub_df.Term.unique())>10:
        ax = gseapy.dotplot(sub_df, title='',cmap='viridis_r', size=20, figsize=(3,6),top_term=30,show_ring=True, ofname=f'../figures/{ds}_{sample_name}_svgi_pw_{x}(supp).png')
    elif len(sub_df.Term.unique())<5:
        ax = gseapy.dotplot(sub_df, title='',cmap='viridis_r', size=30, figsize=(3,8),top_term=30,show_ring=True, ofname=f'../figures/{ds}_{sample_name}_svgi_pw_{x}_large(supp).png')
        ax = gseapy.dotplot(sub_df, title='',cmap='viridis_r', size=100, figsize=(3,1),top_term=30,show_ring=True, ofname=f'../figures/{ds}_{sample_name}_svgi_pw_{x}_small(supp).png')
    else:
        ax = gseapy.dotplot(sub_df, title='',cmap='viridis_r', size=10, figsize=(3,8),top_term=30,show_ring=True, ofname=f'../figures/{ds}_{sample_name}_svgi_pw_{x}_large(supp).png')
        ax = gseapy.dotplot(sub_df, title='',cmap='viridis_r', size=50, figsize=(3,2),top_term=30,show_ring=True, ofname=f'../figures/{ds}_{sample_name}_svgi_pw_{x}_small(supp).png')


op.vis.pattern_LRI(idata,show_SVI=10, spot_size=1)
plt.tight_layout()
plt.savefig(f'../figures/{ds}_{sample_name}_patterns(supp).png', dpi=600,bbox_inches='tight')
plt.close()

op.svi.eva_pattern(idata)
plt.savefig(f'../figures/{ds}_{sample_name}_pattern_metric(supp).png', dpi=600,bbox_inches='tight')
plt.close()



idata = idata[:, idata.var['is_svi']==1]


adata, adata_lri, adata_pattern = op.cell_transform(idata, adata, label=adata.uns['cluster_key'])


with plt.rc_context():
    sc.pl.rank_genes_groups_dotplot(adata_lri, standard_scale='var', show=False)
    plt.savefig(f'../figures/{ds}_{sample_name}_label_SVI(supp).png', bbox_inches="tight")
    plt.close()


with plt.rc_context():
    sc.pl.rank_genes_groups_dotplot(adata_lri, standard_scale='var', n_genes=1, show=False)
    plt.savefig(f'../figures/{ds}_{sample_name}_label_SVI(supp).png', bbox_inches="tight")
    plt.close()
    


with plt.rc_context(): 
    sc.pl.rank_genes_groups_dotplot(adata_pattern, standard_scale='var', n_genes=1, show=False)
    plt.savefig(f'../figures/{ds}_{sample_name}_label_pattern(supp).png', bbox_inches="tight")
    plt.close()



from scipy import stats
pds = []
for i in adata_lri.var_names:
    pds.append(pd.get_dummies(adata.obs[adata.uns['cluster_key']]).corrwith(adata_lri.to_df()[i].astype('float'), method=stats.pointbiserialr)[:1])
pds_df = pd.concat(pds)
pds_df.index = adata_lri.var_names
useful_df = []

df_plot = pd.concat([idata.uns['cell_meta'][['x_global_affine', 'y_global_affine']], pd.get_dummies(idata.uns['cell_meta']['celltype_mapped_refined'])], axis=1)
for i in adata.obs['celltype_mapped_refined'].unique():
    df_sub = pds_df[pds_df[i]>0.5][i].sort_values(ascending=False)
    if len(df_sub) > 0:
        useful_df.append(df_sub)
        markers = df_sub.index.to_numpy()
        corrs = df_sub.values
        cluster_name = df_sub.name
        plt.figure(figsize=(4*(len(markers)+1), 4))
        base = 1
        plt.subplot(1, len(markers)+1, base)

        plt.scatter(df_plot['x_global_affine'],df_plot['y_global_affine'], c=df_plot[cluster_name], s=0.5, cmap='plasma')
        plt.axis('equal')
        plt.axis('off')
        plt.title(f'{cluster_name}\n({sample_name})')
        base += 1

        marker_score = adata_lri.to_df()[markers]
        marker_score = pd.concat([adata.obs[['x_global_affine', 'y_global_affine']], marker_score], axis=1)
        for i in range(len(markers)):
            plt.subplot(1, len(markers)+1, base)
            im=plt.scatter(marker_score['x_global_affine'],marker_score['y_global_affine'], c=marker_score[markers[i]], s=0.5, cmap='plasma')
            plt.colorbar(im,fraction=0.046, pad=0.04)
            plt.axis('equal')
            plt.axis('off')
            plt.title(f'{markers[i]}\ncorr={"%.3f" % corrs[i]}')
            base += 1
        plt.savefig(f'../figures/{ds}_{sample_name}_{cluster_name.replace("/", "-")}(main3J-Morsupp).png', bbox_inches="tight", dpi=300)
        plt.close()


