import matplotlib.pyplot as plt
import seaborn as sns
from spider import SPIDER
op = SPIDER()
import anndata
import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import umap

import os
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
TF_ENABLE_ONEDNN_OPTS = 0
from importlib import resources

sample_name = '151672'
ds = 'DFPLC'
out_f = f'../input_datasets/{ds}/{sample_name}/'
R_path = 'source /etc/profile;module load GCC/11.2.0 OpenMPI/4.1.1 R/4.2.0 Anaconda3/2022.05 R-bundle-Bioconductor/3.15-R-4.2.0;R'
adata = anndata.read_h5ad(f'{out_f}/adata.h5ad')



no_spatalk = False
if len(adata) > 10000:
    no_spatalk=True
idata = op.prep(adata, out_f, R_path, cluster_key=adata.uns['cluster_key'], is_human=adata.uns['is_human'], coord_type=adata.uns['coord_type'], no_spatalk=no_spatalk)
k = idata.uns['cell_meta'][adata.uns['cluster_key']].nunique()
i = 7
j = 50
interface_id = np.where(idata.obs_names == idata.obs[idata.obs['label']=='WM_WM'].index[0])[0]
count_f = f'{out_f}/idata_full_count.csv'
meta_f = f'{out_f}/idata_full_meta.csv'
idata.to_df().to_csv(count_f)
idata.obs[['row', 'col']].to_csv(meta_f)
with resources.path("spider.R_script", "run_spatialPCA.R") as pw_fn:
    os.system(str(f'/bin/bash -c "{R_path} -f /{pw_fn} {meta_f} {count_f} {k} {j} {out_f} {interface_id} {i}"'))


features = pd.read_csv(f'{out_f}/interface_SpatialPCs.csv', index_col=0).T
idata.obsm['spatialPC'] = features.loc[idata.obs_names].to_numpy()
label_df = pd.read_csv(f'{out_f}/refined_interface_label.csv', index_col=0)
idata.obs['label_spatialPCA'] = label_df.loc[idata.obs_names]['clusterlabel_refine'].to_numpy()





umap_3 = umap.UMAP(random_state=52, n_components=3).fit_transform(idata.obsm['spatialPC'])
normalized_df=pd.DataFrame(umap_3)
normalized_df = (normalized_df-normalized_df.min())/(normalized_df.max()-normalized_df.min())
normalized_df.columns=['r', 'g', 'b']
normalized_df[['x', 'y']] = idata.obsm['spatial']
plt.scatter(normalized_df.x, normalized_df.y, c = normalized_df[['r', 'g', 'b']].to_numpy(), s=3)
plt.axis('equal')
plt.axis('off')
plt.savefig(f'../figures/{ds}_{sample_name}_spatialpca_features_rgb(main2B).png', dpi=600,bbox_inches='tight')
plt.close()


g=sns.scatterplot(data=label_df, x='row', y='col', hue = 'clusterlabel_refine', s=6, palette='pastel', linewidth=0)
g.legend(loc='upper left', bbox_to_anchor=(-0.4, 1), ncol=1,frameon =False)
plt.axis('equal')
plt.axis('off')
plt.savefig(f'../figures/{ds}_{sample_name}_spatialpca_cluster(main2C).png', dpi=600,bbox_inches='tight')
plt.close()



dfs = []
for i in idata.obs['label_spatialPCA'].unique():
    dfs.append(idata.uns['cell_meta'].loc[np.unique(idata.obs[idata.obs['label_spatialPCA']==i][['A', 'B']].to_numpy().flatten()), 'layer_guess'].value_counts())
dfs_merged = pd.concat(dfs, axis=1)
dfs_merged.columns=idata.obs['label_spatialPCA'].unique()
dfs_merged = (dfs_merged / dfs_merged.sum()).T
order = np.dot(dfs_merged, (range(1, len(dfs_merged.columns)+1))).argsort()
rename_dict = {}
for i in range(len(dfs_merged.index)):
    rename_dict[dfs_merged.index[order[i]]] = i+1
dfs_merged = dfs_merged.reset_index().rename(columns={'index':'interface cluster'})
dfs_merged['interface cluster'] = dfs_merged['interface cluster'].astype('category').cat.rename_categories(rename_dict).cat.reorder_categories(range(1, len(dfs_merged)+1))
label_df['clusterlabel_refine'] = label_df['clusterlabel_refine'].astype('category').cat.rename_categories(rename_dict).cat.reorder_categories(range(1, len(dfs_merged)+1))
g = dfs_merged.sort_values('interface cluster').plot.bar(x='interface cluster', stacked=True, legend=False)
g.spines[['right', 'top']].set_visible(False)
plt.savefig(f'../figures/{ds}_{sample_name}_spatialpca_cluster_barplot(main2D).png', dpi=600,bbox_inches='tight')
plt.close()


idata.obsm['smooth_pattern_score']  = idata.obsm['spatialPC']
op.traj.paga_default(idata, label='label')
sc.set_figure_params(dpi=300)
idata.obs['label'] = idata.obs['label'].astype('category').cat.rename_categories({
    'WM_WM': 'WM',    
    'Layer4_Layer4': 'Layer4',    
    'Layer1_Layer1': 'Layer1',    
    'Layer5_Layer5': 'Layer5',    
    'Layer6_Layer6': 'Layer6',    
    'Layer2_Layer2': 'Layer2',    
    'Layer3_Layer3': 'Layer3',    
})
with plt.rc_context():
    sc.pl.paga_compare(idata, legend_fontsize=5, frameon=False, size=20, fontweight='normal', legend_loc='on data',
                    title=sample_name, legend_fontoutline=2, show=False)
    plt.savefig(f'../figures/{ds}_{sample_name}_paga_compare(Main2F).png', bbox_inches="tight")
    


idata.obs['label_spatialPCA'] = idata.obs['label_spatialPCA'].astype(str).astype('category')
sc.tl.rank_genes_groups(idata, groupby='label_spatialPCA', method='wilcoxon')
marker_df = op.util.get_marker_df(idata, logfoldchanges_threhold=5)


with plt.rc_context():
    sc.pl.rank_genes_groups_dotplot(idata, n_genes=5, standard_scale='var',)
    # sc.pl.rank_genes_groups_dotplot(idata, n_genes=5,values_to_plot='logfoldchanges', min_logfoldchange=5, vmax=70, vmin=-70, cmap='bwr')
    plt.savefig(f'../figures/{ds}_{sample_name}_spatialpca_marker_dot(supp).png', bbox_inches="tight")


with plt.rc_context():
    sc.pl.spatial(idata, color=list(idata.uns['rank_genes_groups']['names'][:1].flatten()[0]), spot_size=1)
    plt.savefig(f'../figures/{ds}_{sample_name}_spatialpca_marker_exp(supp).png', bbox_inches="tight")


custom, background = op.er.pathway_prep(idata)
marker_df['g'] = marker_df['names']
marker_df['membership'] = 1
marker_df.to_csv(f'../tables/{ds}_{sample_name}_spatialpca_marker.csv')
merged_df_edge, arr = op.er.enrichment(custom, background, marker_df.set_index('names'), groupby='cluster')
merged_df_edge[merged_df_edge["Adjusted P-value"]<=0.05].to_csv(f'../tables/{ds}_{sample_name}_spatialpca_marker_edge_enrichment.csv')
merged_df_edge.to_csv(f'../tables/{ds}_{sample_name}_spatialpca_marker_edge_enrichment_full.csv')
merged_df, arr = op.er.enrichment_interacrtion_gene_df(marker_df.set_index('names'), groupby='cluster')
merged_df[merged_df["Adjusted P-value"]<=0.05].to_csv(f'../tables/{ds}_{sample_name}_spatialpca_marker_node_enrichment.csv')


