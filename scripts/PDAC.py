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

sample_name = 'PDAC_A'
ds = 'PDAC'
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
plt.savefig(f'../figures/{ds}_{sample_name}_metric(main3B).png', dpi=300,bbox_inches='tight')


op.util.adata_moranI(adata, out_f=out_f)
plt.rcParams['font.size'] = 13
merged_df,lri_pw_list,gene_lr_list,gene_pw_list = op.vis.svg_svi_relation(adata, idata, title='PDAC', is_human=adata.uns['is_human'], top=50)
plt.savefig(f'../figures/{ds}_{sample_name}_relation(main4A).png', dpi=600,bbox_inches='tight')





import gseapy
membership=pd.get_dummies(merged_df.set_index('Term')['group']).groupby('Term').sum().astype(str).agg('-'.join, axis=1).reset_index()
for x in membership[0].unique():
    sub_df = merged_df[merged_df.Term.isin(membership[membership[0]==x].Term)]
    if len(sub_df.Term.unique())>20:
        ax = gseapy.dotplot(sub_df, title='',cmap='viridis_r', size=20, figsize=(3,8),top_term=30,show_ring=True, ofname=f'../figures/{ds}_{sample_name}_svgi_pw_{x}.png')
    elif len(sub_df.Term.unique())>15:
        ax = gseapy.dotplot(sub_df, title='',cmap='viridis_r', size=20, figsize=(3,6),top_term=30,show_ring=True, ofname=f'../figures/{ds}_{sample_name}_svgi_pw_{x}.png')
    elif len(sub_df.Term.unique())<5:
        ax = gseapy.dotplot(sub_df, title='',cmap='viridis_r', size=30, figsize=(3,8),top_term=30,show_ring=True, ofname=f'../figures/{ds}_{sample_name}_svgi_pw_{x}_large.png')
        ax = gseapy.dotplot(sub_df, title='',cmap='viridis_r', size=100, figsize=(3,1),top_term=30,show_ring=True, ofname=f'../figures/{ds}_{sample_name}_svgi_pw_{x}_small.png')
    else:
        ax = gseapy.dotplot(sub_df, title='',cmap='viridis_r', size=10, figsize=(3,8),top_term=30,show_ring=True, ofname=f'../figures/{ds}_{sample_name}_svgi_pw_{x}_large.png')
        ax = gseapy.dotplot(sub_df, title='',cmap='viridis_r', size=20, figsize=(3,2),top_term=30,show_ring=True, ofname=f'../figures/{ds}_{sample_name}_svgi_pw_{x}_small.png')


op.vis.pattern_LRI(idata,show_SVI=10, spot_size=10)
plt.tight_layout()
plt.savefig(f'../figures/{ds}_{sample_name}_patterns(supp).png', dpi=600,bbox_inches='tight')


op.svi.eva_pattern(idata)
plt.savefig(f'../figures/{ds}_{sample_name}_pattern_metric(supp).png', dpi=600,bbox_inches='tight')



idata = idata[:, idata.var['is_svi']==1]


adata, adata_lri, adata_pattern = op.cell_transform(idata, adata, label=adata.uns['cluster_key'])


with plt.rc_context():
    sc.pl.rank_genes_groups_dotplot(adata_lri, standard_scale='var', show=False)
    plt.savefig(f'../figures/{ds}_{sample_name}_label_SVI(supp).png', bbox_inches="tight")


with plt.rc_context(): 
    sc.pl.rank_genes_groups_dotplot(adata_pattern, standard_scale='var', n_genes=1, show=False)
    plt.savefig(f'../figures/{ds}_{sample_name}_label_pattern(supp).png', bbox_inches="tight")


st_decon = adata.obsm['deconvolution']
st_decon_org = st_decon.copy()
decon_corr = pd.concat([st_decon, adata_lri.to_df()], axis=1).corr()[adata_lri.var_names].loc[st_decon.columns].T
decon_corr_sup = decon_corr.loc[decon_corr.idxmax()] 
decon_corr = decon_corr[decon_corr_sup.columns[np.argsort(-decon_corr_sup.to_numpy().diagonal())]]
decon_corr_sup = decon_corr.loc[decon_corr.idxmax()] 
plt.rcParams.update({'font.size': 18})
st_decon = st_decon_org.copy()
celtypes = decon_corr_sup.columns
st_decon = pd.concat([idata.uns['cell_meta'][['x', 'y']], st_decon], axis=1)
plt.figure(figsize=(len(celtypes)*4, 8))
base = 1
for i in celtypes:
    plt.subplot(2, len(celtypes), base)
    # sns.scatterplot(st_decon, x='x', y='y', hue=i, s=10, linewidth=0)
    im=plt.scatter(st_decon['x'],st_decon['y'], c=st_decon[i], s=50, cmap='plasma')
    plt.colorbar(im,fraction=0.046, pad=0.04)

    plt.axis('equal')
    plt.axis('off')
    plt.title(i.replace('_', ' '))
    base += 1
    
corrs = decon_corr_sup.to_numpy().diagonal()

marker_score = adata_lri.to_df()[decon_corr_sup.index]
marker_score = pd.concat([adata_lri.obs[['x', 'y']], marker_score], axis=1)
for i in range(len(decon_corr_sup.index)):
    plt.subplot(2, len(celtypes), base)
    im=plt.scatter(marker_score['x'],marker_score['y'], c=marker_score[decon_corr_sup.index[i]].to_numpy(), s=50, cmap='plasma')
    # sns.scatterplot(marker_score, x='x', y='y', hue=i, s=10, linewidth=0)
    plt.colorbar(im,fraction=0.046, pad=0.04)
    plt.axis('equal')
    plt.axis('off')
    plt.title(f'{decon_corr_sup.index[i]} \n corr={"%.3f" % corrs[i]}')
    base += 1

plt.tight_layout()
plt.savefig(f'../figures/{ds}_{sample_name}_decon_SVI(main3C).png', dpi=600,bbox_inches='tight')



