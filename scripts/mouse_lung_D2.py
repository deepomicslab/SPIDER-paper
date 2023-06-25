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


ds = 'mouse_lung'
sample_name = 'D2'
out_f = f'../SPIDER-paper/input_datasets/{ds}/{sample_name}/'
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
plt.savefig(f'../figures/{ds}_{sample_name}_metric(main3partialD).png', dpi=300,bbox_inches='tight')
plt.close()

op.util.adata_moranI(adata, out_f=out_f)
plt.rcParams['font.size'] = 13

merged_df,lri_pw_list,gene_lr_list,gene_pw_list = op.vis.svg_svi_relation(adata, idata, title=sample_name, is_human=adata.uns['is_human'], top=50)
plt.savefig(f'../figures/{ds}_{sample_name}_relation(main4C).png', dpi=600,bbox_inches='tight')
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
    plt.savefig(f'../figures/{ds}_{sample_name}_label_SVI1(supp).png', bbox_inches="tight")
    plt.close()


with plt.rc_context(): 
    sc.pl.rank_genes_groups_dotplot(adata_pattern, standard_scale='var', n_genes=1, show=False)
    plt.savefig(f'../figures/{ds}_{sample_name}_label_pattern(supp).png', bbox_inches="tight")
    plt.close()



from scipy import stats
pds = []
for i in adata_lri.var_names:
    pds.append(pd.get_dummies(adata.obs['annotation']).corrwith(adata_lri.to_df()[i].astype('float'), method=stats.pointbiserialr)[:1])
pds_df = pd.concat(pds)
pds_df.index = adata_lri.var_names
adata.obs[['row', 'col']] = adata.obs[['row', 'col']].astype(int)
adata_lri.obs[['row', 'col']] = adata_lri.obs[['row', 'col']].astype(int)
useful_df = []
plt.rcParams.update({'font.size': 14})
for cluster_i in adata.obs['annotation'].unique():
    df_sub = pds_df[pds_df[cluster_i]>0.5][cluster_i].sort_values(ascending=False)
    if len(df_sub) > 0:
        useful_df.append(df_sub)
        markers = df_sub.index.to_numpy()[:10]
        corrs = df_sub.values[:10]
        cluster_name = df_sub.name
        plt.figure(figsize=(4*(len(markers)+1), 4))
        base = 1
        plt.subplot(1, len(markers)+1, base)

        plt.scatter(adata.obs['row'],adata.obs['col'], c=adata.obs['annotation']==cluster_i, s=0.5, cmap='plasma')
        plt.axis('equal')
        plt.axis('off')
        plt.title(f'{cluster_name}\n({sample_name})')
        base += 1

        marker_score = adata_lri.to_df()[markers]
        marker_score = pd.concat([adata_lri.obs[['row', 'col']], marker_score], axis=1)
        for i in range(len(markers)):
            plt.subplot(1, len(markers)+1, base)
            im=plt.scatter(marker_score['row'],marker_score['col'], c=marker_score[markers[i]], s=0.5, cmap='plasma')
            plt.colorbar(im,fraction=0.046, pad=0.04)
            plt.axis('equal')
            plt.axis('off')
            plt.title(f'{markers[i]}\ncorr={"%.3f" % corrs[i]}')
            base += 1
        plt.savefig(f'../figures/{ds}_{sample_name}_{cluster_name}(main3G).png', dpi=300,bbox_inches='tight')
        plt.close()


arr = []
for cluster_i in adata.obs['annotation'].unique():
    df_sub = pds_df[pds_df[cluster_i]>0.5][cluster_i].sort_values(ascending=False)
    for j in df_sub.index:
        arr.append([cluster_i, j])
celltype_lris = pd.DataFrame(arr, columns=['celltype', 'lri']).set_index('lri')
merged_df, arr = op.er.enrichment_interacrtion_gene_df(celltype_lris, groupby='celltype',  is_human=False)
merged_df = merged_df[merged_df['Adjusted P-value'] < 0.05]
merged_df[merged_df.Term.str.contains('adhesion')]


adata.obsm['spatial'] = adata.obsm['spatial'].astype(int)
sc.pl.spatial(adata, color=['Scgb1a1', 'Scgb3a2', 'Cyp2f2', 'Hp'], spot_size=1.5)


corr_df = pd.concat([adata.to_df()[['Scgb1a1', 'Scgb3a2', 'Cyp2f2', 'Hp']].loc[adata_lri.obs_names], adata_lri.to_df()], axis=1).corr(method='pearson')
corr_df_asyn =  corr_df.loc[adata_lri.var_names, ['Scgb1a1', 'Scgb3a2', 'Cyp2f2', 'Hp']]
plt.rcParams.update({'font.size': 15})

corrs = corr_df_asyn.mean(1).sort_values(ascending=False).values[:10]
lris = corr_df_asyn.mean(1).sort_values(ascending=False).index[:10]
marker_score = adata_lri.to_df()[lris]
adata.obs[['row', 'col']] = adata.obs[['row', 'col']].astype(int)
adata_lri.obs[['row', 'col']] = adata_lri.obs[['row', 'col']].astype(int)

plt.figure(figsize=(44, 4))
base=1
plt.subplot(1, 11, base)
im=plt.scatter(adata.obs['row'],adata.obs['col'], c=adata.to_df()['Scgb1a1'], s=0.05, cmap='plasma',vmin=0)
plt.colorbar(im,fraction=0.046, pad=0.04)
plt.axis('equal')
plt.axis('off')
plt.title('Bronchioli\nScgb1a1')
base+=1
for i in range(len(lris)):
    plt.subplot(1, 11, base)
    im=plt.scatter(adata_lri.obs['row'],adata_lri.obs['col'], c=marker_score[lris[i]], s=0.05, cmap='plasma')
    plt.colorbar(im,fraction=0.046, pad=0.04)
    plt.axis('equal')
    plt.axis('off')
    plt.title(f'{lris[i]} \n corr={"%.3f" % corrs[i]}')
    base += 1
plt.tight_layout()
plt.savefig(f'../figures/{ds}_{sample_name}_BronchioliMarker(main3F).png', dpi=600,bbox_inches='tight')
plt.close()


decon = pd.read_csv(f'{out_f}/decon.csv', index_col=0)
corr_df = pd.concat([decon.T, adata_lri.to_df().loc[decon.columns]], axis=1).corr(method='pearson')
corr_df_asyn =  corr_df.loc[adata_lri.var_names, decon.index]
mask = corr_df_asyn.loc[corr_df_asyn.idxmax()].to_numpy().diagonal() > 0.3
decon_meta = pd.concat([decon.T, adata_lri.obs.loc[decon.columns]], axis=1)
decon_meta[['row', 'col']] = decon_meta[['row', 'col']].astype(int)


corr_df_asyn.loc[corr_df_asyn.idxmax()[mask], corr_df_asyn.idxmax()[mask].index]


plt.rcParams.update({'font.size': 15})

corrs = corr_df_asyn[corr_df_asyn['Adventitial Fibroblast']>0.3]['Adventitial Fibroblast'].sort_values(ascending=False).values[:10]
lris = corr_df_asyn[corr_df_asyn['Adventitial Fibroblast']>0.3]['Adventitial Fibroblast'].sort_values(ascending=False).index[:10]
marker_score = adata_lri.to_df().loc[decon.columns, lris]

plt.figure(figsize=(44, 4))
base=1
plt.subplot(1, 11, base)
cell_meta_sub=decon_meta[decon_meta['Adventitial Fibroblast']>0].sort_values('Adventitial Fibroblast')
im=plt.scatter(decon_meta['row'],decon_meta['col'], c=decon_meta['Adventitial Fibroblast'], s=1, cmap='plasma',vmin=0, edgecolors='none')
plt.scatter(cell_meta_sub['row'],cell_meta_sub['col'], c=cell_meta_sub['Adventitial Fibroblast'], s=3, cmap='plasma',vmin=0, edgecolors='none')
plt.colorbar(im,fraction=0.046, pad=0.04)
plt.axis('equal')
plt.axis('off')
plt.title(f'Adventitial\nFibroblast')
base+=1
for i in range(len(lris)):
    plt.subplot(1, 11, base)
    im=plt.scatter(decon_meta['row'],decon_meta['col'], c=marker_score[lris[i]], s=1, cmap='plasma', edgecolors='none')
    plt.colorbar(im,fraction=0.046, pad=0.04)
    plt.axis('equal')
    plt.axis('off')
    plt.title(f'{lris[i]} \n corr={"%.3f" % corrs[i]}')
    base += 1
plt.tight_layout()
plt.savefig(f'../figures/{ds}_{sample_name}_AdventitialFibroblast(main3H).png', dpi=600,bbox_inches='tight')
plt.close()


# plt.rcParams.update({'font.size': 22})

corrs = corr_df_asyn[corr_df_asyn['Alveolar Fibroblast']>0.3]['Alveolar Fibroblast'].sort_values(ascending=False).values[:10]
lris = corr_df_asyn[corr_df_asyn['Alveolar Fibroblast']>0.3]['Alveolar Fibroblast'].sort_values(ascending=False).index[:10]
marker_score = adata_lri.to_df().loc[decon.columns, lris]

plt.figure(figsize=(44, 4))
base=1
plt.subplot(1, 11, base)
im=plt.scatter(decon_meta['row'],decon_meta['col'], c=decon_meta['Alveolar Fibroblast'], s=0.05, cmap='plasma',vmin=0)
plt.colorbar(im,fraction=0.046, pad=0.04)
plt.axis('equal')
plt.axis('off')
plt.title(f'Alveolar\nFibroblast')
base+=1
for i in range(len(lris)):
    plt.subplot(1, 11, base)
    im=plt.scatter(decon_meta['row'],decon_meta['col'], c=marker_score[lris[i]], s=0.05, cmap='plasma')
    plt.colorbar(im,fraction=0.046, pad=0.04)
    plt.axis('equal')
    plt.axis('off')
    plt.title(f'{lris[i]} \n corr={"%.3f" % corrs[i]}')
    base += 1
plt.tight_layout()
plt.savefig(f'../figures/{ds}_{sample_name}_AlveolarFibroblast(supp).png', dpi=600,bbox_inches='tight')
plt.close()


plt.rcParams.update({'font.size': 22})

corrs = corr_df_asyn[corr_df_asyn['Ciliated']>0.3]['Ciliated'].sort_values(ascending=False).values[:10]
lris = corr_df_asyn[corr_df_asyn['Ciliated']>0.3]['Ciliated'].sort_values(ascending=False).index[:10]
marker_score = adata_lri.to_df().loc[decon.columns, lris]

plt.figure(figsize=(44, 4))
base=1
plt.subplot(1, 11, base)
im=plt.scatter(decon_meta['row'],decon_meta['col'], c=decon_meta['Ciliated'], s=0.05, cmap='plasma',vmin=0)
plt.colorbar(im,fraction=0.046, pad=0.04)
plt.axis('equal')
plt.axis('off')
plt.title(f'Ciliated')
base+=1
for i in range(len(lris)):
    plt.subplot(1, 11, base)
    im=plt.scatter(decon_meta['row'],decon_meta['col'], c=marker_score[lris[i]], s=0.05, cmap='plasma')
    plt.colorbar(im,fraction=0.046, pad=0.04)
    plt.axis('equal')
    plt.axis('off')
    plt.title(f'{lris[i]} \n corr={"%.3f" % corrs[i]}')
    base += 1
plt.tight_layout()
plt.savefig(f'../figures/{ds}_{sample_name}_Ciliated(supp).png', dpi=600,bbox_inches='tight')
plt.close()


arr = []
for cluster_i in corr_df_asyn.columns:
    df_sub = corr_df_asyn[corr_df_asyn[cluster_i]>0.3][cluster_i].sort_values(ascending=False)
    for j in df_sub.index:
        arr.append([cluster_i, j])
celltype_lris = pd.DataFrame(arr, columns=['celltype', 'lri']).set_index('lri')
merged_df, arr = op.er.enrichment_interacrtion_gene_df(celltype_lris, groupby='celltype',  is_human=False)
merged_df[merged_df['celltype']=='Adventitial Fibroblast']


