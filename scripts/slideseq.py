#!/usr/bin/env /home/lishiying/.conda/envs/spider1/bin/python

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

sample_name = 'mousebrain'
ds = 'slide-seq-v2'
out_f = f'../input_datasets/{ds}/{sample_name}/'
adata = anndata.read_h5ad(f'{out_f}/adata.h5ad')

R_path = 'source /etc/profile;module load GCC/11.2.0 OpenMPI/4.1.1 R/4.2.0 Anaconda3/2022.05 R-bundle-Bioconductor/3.15-R-4.2.0;R'


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
merged_df,lri_pw_list,gene_lr_list,gene_pw_list = op.vis.svg_svi_relation(adata, idata, title='mouse brain', is_human=adata.uns['is_human'], top=50)
plt.savefig(f'../figures/{ds}_{sample_name}_relation(main4B).png', dpi=600,bbox_inches='tight')
plt.close()


import gseapy
membership=pd.get_dummies(merged_df.set_index('Term')['group']).groupby('Term').sum().astype(str).agg('-'.join, axis=1).reset_index()
for x in membership[0].unique():
    sub_df = merged_df[merged_df.Term.isin(membership[membership[0]==x].Term)]
    if len(sub_df.Term.unique())>20:
        ax = gseapy.dotplot(sub_df, title='',cmap='viridis_r', size=20, figsize=(3,8),top_term=30,show_ring=True, ofname=f'../figures/{ds}_{sample_name}_svgi_pw_{x}(supp).png')
    elif len(sub_df.Term.unique())>15:
        ax = gseapy.dotplot(sub_df, title='',cmap='viridis_r', size=20, figsize=(3,6),top_term=30,show_ring=True, ofname=f'../figures/{ds}_{sample_name}_svgi_pw_{x}(supp).png')
    elif len(sub_df.Term.unique())<5:
        ax = gseapy.dotplot(sub_df, title='',cmap='viridis_r', size=30, figsize=(3,8),top_term=30,show_ring=True, ofname=f'../figures/{ds}_{sample_name}_svgi_pw_{x}_large(supp).png')
        ax = gseapy.dotplot(sub_df, title='',cmap='viridis_r', size=100, figsize=(3,1),top_term=30,show_ring=True, ofname=f'../figures/{ds}_{sample_name}_svgi_pw_{x}_small(supp).png')
    else:
        ax = gseapy.dotplot(sub_df, title='',cmap='viridis_r', size=10, figsize=(3,8),top_term=30,show_ring=True, ofname=f'../figures/{ds}_{sample_name}_svgi_pw_{x}_large(supp).png')
        ax = gseapy.dotplot(sub_df, title='',cmap='viridis_r', size=20, figsize=(3,2),top_term=30,show_ring=True, ofname=f'../figures/{ds}_{sample_name}_svgi_pw_{x}_small(supp).png')


op.vis.pattern_LRI(idata,show_SVI=10, spot_size=1)
plt.tight_layout()
plt.savefig(f'../figures/{ds}_{sample_name}_patterns(supp).png', dpi=600,bbox_inches='tight')
plt.close()


op.svi.eva_pattern(idata)
plt.savefig(f'../figures/{ds}_{sample_name}_pattern_metric(main3P).png', dpi=600,bbox_inches='tight')
plt.close()



idata = idata[:, idata.var['is_svi']==1]


adata, adata_lri, adata_pattern = op.cell_transform(idata, adata, label=adata.uns['cluster_key'])


with plt.rc_context():
    sc.pl.rank_genes_groups_dotplot(adata_lri, standard_scale='var', show=False)
    plt.savefig(f'../figures/{ds}_{sample_name}_label_SVI(supp).png', bbox_inches="tight")
    plt.close()


with plt.rc_context(): 
    sc.pl.rank_genes_groups_dotplot(adata_pattern, standard_scale='var', n_genes=1, show=False)
    plt.savefig(f'../figures/{ds}_{sample_name}_label_pattern(supp).png', bbox_inches="tight")
    plt.close()



region = pd.get_dummies(idata.uns['cell_meta'].loc[adata.obs_names]['cluster']).loc[adata_pattern.obs_names]
from scipy import stats
pds = []
for i in adata_pattern.var_names:
    pds.append(region.corrwith(adata_pattern.to_df()[i].astype('float'), method=stats.pointbiserialr)[:1])
pds_df = pd.concat(pds)
pds_df.index = adata_pattern.var_names
useful_df = []
for cluster_i in adata.obs['cluster'].unique():
    df_sub = pds_df[pds_df[cluster_i]>0.3][cluster_i].sort_values(ascending=False)
    if len(df_sub) > 0:
        useful_df.append(df_sub)

patterns = []
cluster = []
cluster_id = []
corrs = []
cluster_id_map = {
    'DentatePyramids': 'Dentate Pyramids',
    'Endothelial_Tip': 'Endothelial Tip',
    'Oligodendrocytes': 'Oligodendrocytes',
    'Astrocytes': 'Astrocytes',
    'Ependymal': 'Ependymal'
}
for x in useful_df:
    patterns.append(int(x.idxmax()))
    corrs.append(x.max())
    cluster_id.append(x.name)
    cluster.append(cluster_id_map[x.name])
    


useful_df = []
adata_lri.obsm['spatial'] = adata[adata_lri.obs_names].obsm['spatial']
adata_pattern.obsm['spatial'] = adata[adata_lri.obs_names].obsm['spatial']

plt.figure(figsize=(5*4, 5*len(patterns)))
base = 1
for i in range(len(patterns)):
    corr =  corrs[i]
    cluster_name = cluster[i]
    cluster_i = cluster_id[i]

    plt.subplot(len(patterns), 4, base)

    plt.scatter(adata.obsm['spatial'][:,0],adata.obsm['spatial'][:,1], c=adata.obs['cluster']==cluster_i, s=3, cmap='plasma', edgecolors='none')
    # plt.scatter(adata[adata.obs['cluster']==cluster_id].obsm['spatial'][:,0],adata[adata.obs['cluster']==cluster_id].obsm['spatial'][:,1], c=adata[adata.obs['cluster']==cluster_name].obs['cluster']==cluster_name, s=5, cmap='plasma', vmin=0, edgecolors='none')
    plt.axis('equal')
    plt.axis('off')
    plt.title(f'{cluster_name}\n({sample_name})')
    base += 1
    
    plt.subplot(len(patterns), 4, base)
    im=plt.scatter(adata_pattern.obsm['spatial'][:,0], adata_pattern.obsm['spatial'][:,1], c=adata_pattern.to_df()[str(patterns[i])], s=3, cmap='plasma', edgecolors='none')
    plt.colorbar(im,fraction=0.046, pad=0.04)
    plt.axis('equal')
    plt.axis('off')
    plt.title(f'Pattern {patterns[i]} - {np.sum(idata.var.label == patterns[i])} SVIs\ncorr={"%.3f" % corr}')
    base += 1
    
    svis = idata.var[idata.var.label == patterns[i]].sort_values(f'pattern_correlation_{patterns[i]}', ascending=False)
    svi_to_plot = svis.index.to_numpy()[:2]
    print(svi_to_plot)
    memberships = svis[f'pattern_correlation_{patterns[i]}'].to_numpy()[:2]
    for j in range(2):
        plt.subplot(len(patterns), 4, base)
        im=plt.scatter(adata_lri.obsm['spatial'][:,0], adata_lri.obsm['spatial'][:,1], c=adata_lri.to_df()[svi_to_plot[j]], s=3, cmap='plasma', edgecolors='none',linewidths=0)
        plt.axis('equal')
        plt.axis('off')
        plt.colorbar(im,fraction=0.046, pad=0.04)
        plt.title(f'{svi_to_plot[j]}\n corr={"%.3f" % memberships[j]}')
        base += 1
plt.savefig(f'../figures/{ds}_{sample_name}_cluster_pattern(main3O).png', dpi=100,bbox_inches='tight')
plt.close()


group = idata.var[idata.var['is_svi']==1]
group['g'] = group.index
group['membership'] = 1
custom, background = op.er.pathway_prep(idata,is_human=False)
merged_df_edge, arr = op.er.enrichment(custom, background, group, groupby='label')



merged_df_edge[merged_df_edge.group=='11']


enrich_df_gene_reactome, arr = op.er.enrichment_interacrtion_gene_df(group, groupby='label', is_human=False, custom_pathwaydb=['Reactome_2022'])
enrich_df_gene_reactome['ratio'] = np.array([int(x.split('/')[0])/int(x.split('/')[1]) for x in enrich_df_gene_reactome.Overlap])
enrich_df_gene_reactome['pattern_rep'] = pd.Categorical(enrich_df_gene_reactome.label).rename_categories({
    '11': 'Dentate Pyramids (11)',
    '7': 'Oligodendrocytes (7)', 
    '0': 'Astrocytes (0)', 
    '5': 'Ependymal (5)', 
    '6': 'Endothelial Tip (6)', 
})


enrich_df_gene_reactome['Term_old'] = enrich_df_gene_reactome['Term']
enrich_df_gene_reactome['Term'] = [s.rsplit(' ', 1)[0] for s in enrich_df_gene_reactome['Term_old']]
enrich_df_gene_reactome['Term'] = enrich_df_gene_reactome['Term'].replace('p130Cas Linkage To MAPK Signaling For Integrins', 'MAPK Signaling For Integrins (1)')
enrich_df_gene_reactome['Term'] = enrich_df_gene_reactome['Term'].replace('GRB2:SOS Provides Linkage To MAPK Signaling For Integrins', 'MAPK signaling For integrins (2)')
enrich_df_gene_reactome['Term'] = enrich_df_gene_reactome['Term'].replace('Regulation Of Commissural Axon Pathfinding By SLIT And ROBO', 'Regulation Of Commissural Axon Pathfinding')
enrich_df_gene_reactome = enrich_df_gene_reactome[enrich_df_gene_reactome['Adjusted P-value']<0.05]
op.vis.enrichment(enrich_df_gene_reactome[(enrich_df_gene_reactome.label.isin(np.array(patterns, dtype=str))) & (enrich_df_gene_reactome['ratio'] > 0.12)], x_key='pattern_rep', size=5, figsize=(3,7), top_term=5, save=f'../figures/{ds}_{sample_name}_reactome_small(mainQ).png')
op.vis.enrichment(enrich_df_gene_reactome[(enrich_df_gene_reactome.label.isin(np.array(patterns, dtype=str))) & (enrich_df_gene_reactome['ratio'] > 0.12)], x_key='pattern_rep', size=5, figsize=(3,10), top_term=5, save=f'../figures/{ds}_{sample_name}_reactome_large(mainQ).png')
op.vis.enrichment(enrich_df_gene_reactome[(enrich_df_gene_reactome.label.isin(np.array(patterns, dtype=str))) & (enrich_df_gene_reactome['ratio'] > 0.12)], x_key='pattern_rep', size=5, figsize=(3,7), top_term=5)





