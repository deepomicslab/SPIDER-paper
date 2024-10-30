import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator
from spider import SPIDER
op=SPIDER()

def compare_sdm(idata, svi_df_strict, dm_expand):
    sdm_unique = dm_expand.loc[(list(set(dm_expand.index) - set(svi_df_strict.index)))].sort_values('fdr').index
    sp_unique = svi_df_strict.loc[(list(set(svi_df_strict.index) - set(dm_expand.index)))].max(axis=1).sort_values().index
    shared = list(set(dm_expand.index) & set(svi_df_strict.index))
    methods = np.array(['moranI', 'gearyC', 'SOMDE', 'nnSVG'])[np.isin(['SOMDE', 'nnSVG', 'gearyC', 'moranI'],list(idata.uns.keys()))]
    print(f'evaluating with {methods}')
    dfs = []
    metrics = []
    for i in methods:
        if i == 'gearyC':
            dfs.append(-idata.uns['gearyC'][['C']])
            metrics.append("Geary\nC (rev.)")
        elif i == 'moranI':
            dfs.append(idata.uns['moranI'][['I']]),
            metrics.append("Moran\nI")
        elif i == 'SOMDE':
            dfs.append(idata.uns['SOMDE'].set_index('g')['FSV']),
            metrics.append("SOMDE\nFSV") 
            dfs.append(idata.uns['SOMDE'].set_index('g')['LLR']),
            metrics.append("SOMDE\nLR") 
        elif i == 'nnSVG':
            dfs.append(idata.uns['nnSVG']['LR_stat']),
            metrics.append("nnSVG\nLR") 
    metrics.append('TF corr')
    pairs = []
    for i in metrics:
        pairs.append( ((i, 'SPIDER'), (i, 'Excluded')) )
        pairs.append( ((i, 'SpatialDM'), (i, 'Excluded')) )
        pairs.append( ((i, 'SPIDER'), (i, 'SpatialDM')) )
    df = pd.concat(dfs, axis=1)

    normalized_df=(df-df.min())/(df.max()-df.min())
    normalized_df = normalized_df.fillna(0)
    df = pd.DataFrame(idata.uns['tf_corr'].max(axis=1), columns = ['TF corr'], index=idata.var_names)
    normalized_df = pd.concat([normalized_df, df])
    normalized_df.columns=metrics

    normalized_df['Category'] = 'Excluded'
    normalized_df.loc[sdm_unique, 'Category'] = 'SpatialDM'
    normalized_df.loc[sp_unique, 'Category'] = 'SPIDER'
    normalized_df.loc[shared, 'Category'] = 'SPIDER'
    normalized_df = pd.concat([normalized_df, normalized_df.loc[shared].replace('SPIDER', 'SpatialDM')], axis=0)

    plt.rcParams.update({'font.size': 18})
    plt.figure(figsize=(9, 6))

    normalized_df = normalized_df.melt(id_vars='Category', value_vars=metrics, var_name='Metric')
    ax =sns.boxplot(data=normalized_df,x='Metric',y='value', hue='Category', 
                    palette={'SpatialDM':'#313695', 'SPIDER': '#D53E4F', 'Excluded': '#A578D8'}, 
                    width=0.8, hue_order=['SPIDER', 'SpatialDM', 'Excluded'])
    ax.legend(loc='lower left',ncol=1, bbox_to_anchor=(1.1, -0.05), frameon=False)

    annot = Annotator(ax, pairs, data=normalized_df, x='Metric',y='value', hue='Category', hue_order=['SPIDER', 'SpatialDM', 'Excluded'])
    annot.configure(test='Mann-Whitney-gt',comparisons_correction="BH", correction_format="replace")
    annot.apply_and_annotate()
    ax.set_ylabel('')    
    ax.set_xlabel('')
    plt.show()
    
def compare_sdm_tf(idata, dm_expand, ds, sample_name):
    corr_df = pd.DataFrame(idata.uns['tf_corr'].max(axis=1), columns = ['TF corr'], index=idata.var_names)
    corr_df['count'] = idata.var['tf_support_count']
    tf_supported = idata.var.query('tf_support_count >= 1').index.tolist()
    corr_df_spider = corr_df.loc[tf_supported]
    corr_df_spider['Category'] = 'SPIDER-TF'
    corr_df_sdm = corr_df.loc[dm_expand.index]
    corr_df_sdm['Category'] = 'SpatialDM'
    plot_corr_df = pd.concat([corr_df_spider, corr_df_sdm]) 
    plot_corr_df['Sample'] = sample_name
    plot_corr_df['Dataset'] = ds
    
    import seaborn as sns
    from statannotations.Annotator import Annotator
    pairs = [((sample_name, 'SPIDER-TF'), (sample_name, 'SpatialDM'))]
    plt.figure(figsize=(3, 6))
    ax =sns.boxplot(data=plot_corr_df,x='Sample',y='count', hue='Category', medianprops=dict(color="lightgrey"),
                    palette={'SpatialDM':'#313695', 'SPIDER-TF': '#F1BFB7',}, 
                    width=0.8, hue_order=['SPIDER-TF', 'SpatialDM'])
    ax.legend(loc='lower left',ncol=1, bbox_to_anchor=(1.1, -0.05), frameon=False)

    annot = Annotator(ax, pairs, data=plot_corr_df, x='Sample',y='count', hue='Category', hue_order=['SPIDER-TF', 'SpatialDM'])
    annot.configure(test='Mann-Whitney-gt',comparisons_correction="BH", correction_format="replace")
    annot.apply_and_annotate()
    ax.set_ylabel('')    
    ax.set_xlabel('')
    plt.title('Number\nof Supporting svTFs')
    
def top_supported_svi(idata):
    corr_df = pd.DataFrame(idata.uns['tf_corr'].max(axis=1), columns = ['TF corr'], index=idata.var_names)
    corr_df['count'] = idata.var['tf_support_count']
    corr_df['SPIDER'] = 0
    svi_df, svi_df_strict = op.svi.combine_SVI(idata,threshold=0.01)
    corr_df.loc[svi_df_strict.index, 'SPIDER'] = 1
    corr_df.replace(-1, 0, inplace=True)
    corr_df_sub = corr_df[corr_df['SPIDER'] == 1].sort_values('count', ascending=False)
    corr_df_sub = corr_df_sub[corr_df_sub['count'] > 0]
    corr_df_sub['count'].plot.bar(color='#D53E4F', figsize=(10, 5))
    plt.show()
    return corr_df_sub

def go_svi_tf(idata, svis):
    top_5 = svis
    tfs_lris = []
    for lri in top_5:
        support_tf = idata.uns['tf_corr'].loc[lri].index[idata.uns['tf_corr'].loc[lri].values>0.3]
        support_tf = [x.replace(lri, '') for x in support_tf]
        tfs_lris.append(support_tf)

    import gseapy
    is_human=True
    organism = 'Human' if is_human else 'Mouse'
    pathway_db = ['KEGG_2021_Human' if is_human else 'KEGG_2019_Mouse']
    dfs = []
    drop_illness=['Cushing syndrome', 'infection', 'Amoebiasis', 'atherosclerosis', 'cardiomyopathy', 'leukemia', 'Leishmaniasis', 'thyroid', 'Viral', 'diabetic', 'diabetes', 'Asthma', 'arthritis', 'Graft']
    for lri_index in range(len(top_5)):
        lri = top_5[lri_index]
        try:
            support_tf = tfs_lris[lri_index]
            genes = support_tf + lri.split('_')
            
            enr_res = gseapy.enrichr(gene_list=genes,
                            organism=organism,
                            gene_sets=pathway_db)
            enr_res.res2d['LRI'] = lri.replace('_', '-')
            dfs.append(enr_res.res2d)
        except Exception as e:
            print(e)
            continue
    merged_df = pd.concat(dfs)
    merged_df = merged_df[merged_df['Adjusted P-value']<0.05]
    if len(drop_illness) != 0:
        merged_df = merged_df[~merged_df.Term.str.contains('|'.join(drop_illness))]
    sub_df = merged_df.sort_values('Adjusted P-value').groupby('LRI').head(5)
    op.vis.enrichment(sub_df, x_key='LRI', size=10, figsize=(3.2, 5))
    sub_df = sub_df[['LRI', 'Term', 'Adjusted P-value', 'Genes']].sort_values(['LRI', 'Adjusted P-value'])
    return sub_df

def tf_venn(idata, dm_expand):
    from matplotlib_venn import venn2,venn2_circles
    from matplotlib import pyplot as plt
    plt.figure(figsize=(4,4))
    setdm = dm_expand.index
    tf_supported = idata.var.query('tf_support_count >= 1').index.tolist()
    settf = tf_supported
    out=venn2(subsets=[set(settf), set(setdm)], set_labels = ('SPIDER-TF', 'SpatialDM'), set_colors=('#DE614D', 'blue'))
    
def average_lr(idata, dm_expand):
    corr_df = pd.DataFrame(idata.uns['tf_corr'].max(axis=1), columns = ['TF corr'], index=idata.var_names)
    corr_df['count'] = idata.var['tf_support_count']
    corr_df['SPIDERTF'] = 0
    tf_supported = idata.var.query('tf_support_count >= 1').index.tolist()
    corr_df.loc[tf_supported, 'SPIDERTF'] = 1
    corr_df['SPIDER'] = 0
    svi_df, svi_df_strict = op.svi.combine_SVI(idata,threshold=0.01)
    corr_df.loc[svi_df_strict.index, 'SPIDER'] = 1
    corr_df['SpatialDM'] = 0
    corr_df.loc[dm_expand.index, 'SpatialDM'] = 1
    corr_df.replace(-1, 0, inplace=True)
    corr_df['ligand'] = corr_df.index.str.split('_').str[0]
    corr_df['receptor'] = corr_df.index.str.split('_').str[1]
    arr_ligand = [corr_df.query('SpatialDM == 1')['receptor'].value_counts().mean(),
        corr_df.query('SPIDER == 1')['receptor'].value_counts().mean(),  
        corr_df.query('SPIDERTF == 1')['receptor'].value_counts().mean(), corr_df['receptor'].value_counts().mean()]
    arr_receptor = [corr_df.query('SpatialDM == 1')['ligand'].value_counts().mean(), 
                    corr_df.query('SPIDER == 1')['ligand'].value_counts().mean(), 
                    corr_df.query('SPIDERTF == 1')['ligand'].value_counts().mean(), corr_df['ligand'].value_counts().mean()]
    df = pd.DataFrame([arr_ligand, arr_receptor], columns=['SpatialDM', 'SPIDER', 'SPIDER-TF', 'all'], index=['receptor', 'ligand']).T
    df_melt = df.reset_index().melt(id_vars='index', var_name='Metric', value_name='value')
    plt.figure(figsize=(5, 6))
    sns.barplot(data=df_melt, x='Metric', y='value', hue='index', hue_order=['SPIDER', 'SpatialDM', 'SPIDER-TF'],
                palette={'SpatialDM':'#313695', 'SPIDER': '#D53E4F', 'SPIDER-TF': '#F1BFB7', 'all': '#A578D8'})
    plt.ylabel('Average\nNumber of\LRIs', fontsize=18)
    plt.xlabel('')
    
def plot_top_tf(idata, lris):
    cols = []
    tfs = []
    top_corrs = []
    for lri in lris:
        _, top_corr = idata.uns['tf_corr'].loc[lri].sort_values(ascending=False)[:1].reset_index().to_numpy()[0]
        top_corrs.append(top_corr)
    # sort lris based on top_corrs
    lris = [x for _, x in sorted(zip(top_corrs, lris), key=lambda pair: pair[0], reverse=True)]

    for lri in lris:
        tf, corr = idata.uns['tf_corr'].loc[lri].sort_values(ascending=False)[:1].reset_index().to_numpy()[0]
        if lri+'-'+tf in idata.obsm['tf_score'].columns:
            tfs.append(lri+'-'+tf)
            cols.append(f'{tf}\ncorr={"%.3f" % corr}')
    for lri in lris:
        tf, corr = idata.uns['tf_corr'].loc[lri].sort_values(ascending=False)[1:2].reset_index().to_numpy()[0]
        if lri+'-'+tf in idata.obsm['tf_score'].columns:
            tfs.append(lri+'-'+tf)
            cols.append(f'{tf}\ncorr={"%.3f" % corr}')
    corr_df_tf = idata.obsm['tf_score'][tfs]
    corr_df_tf = (corr_df_tf - corr_df_tf.min()) / (corr_df_tf.max() - corr_df_tf.min())
    corr_df_tf.columns = cols
    # drop duplicated columns
    corr_df_tf = corr_df_tf.loc[:,~corr_df_tf.columns.duplicated()]
    corr_df_tf['x'] = idata.obsm['spatial'][:, 0]
    corr_df_tf['y'] = idata.obsm['spatial'][:, 1]
    lri_df = idata.to_df()[lris]
    lri_df[['x', 'y']] = idata.obsm['spatial']

    import matplotlib.pyplot as plt 
    plt.figure(figsize=(len(lris)*2, 6))
    plt.rcParams['font.size'] = 11
    base = 1
    for i in lris:
        plt.subplot(3, len(lris), base)
        # sns.scatterplot(st_decon, x='x', y='y', hue=i, s=10, linewidth=0)
        lri_df = lri_df.sort_values(i)
        im=plt.scatter(lri_df['x'],lri_df['y'], c=lri_df[i], s=1, cmap='plasma')
        plt.colorbar(im,fraction=0.046, pad=0.04)
        plt.axis('equal')
        plt.axis('off')
        plt.title(i)
        base += 1
        
    for i in cols:
        plt.subplot(3, len(lris), base)
        corr_df_tf = corr_df_tf.sort_values(i)
        im=plt.scatter(corr_df_tf['x'],corr_df_tf['y'], c=corr_df_tf[i], s=1, cmap='plasma')
        # sns.scatterplot(marker_score, x='x', y='y', hue=i, s=10, linewidth=0)
        plt.colorbar(im,fraction=0.046, pad=0.04)
        plt.axis('equal')
        plt.axis('off')
        plt.title(i)
        base += 1
    plt.subplots_adjust(wspace=0.5, hspace=0.4) 
    plt.show()
    
def celltype_correlation(adata_svi, celltype_dummy):
    decon_corr = pd.concat([celltype_dummy, adata_svi.to_df()], axis=1).corr()[adata_svi.var_names].loc[celltype_dummy.columns].T
    decon_corr_melt = decon_corr.reset_index().melt(id_vars='index')
    decon_corr_melt.columns = ['SVI', 'celltype', 'correlation']
    decon_corr_melt_sub = decon_corr_melt.groupby('celltype').apply(lambda x: x.nlargest(5, 'correlation')).reset_index(drop=True)
    return decon_corr_melt, decon_corr_melt_sub

def plot_correlation(spider_corr_sub, sdm_corr_sub):
    sdm_corr_sub['method'] = 'SpatialDM'
    spider_corr_sub['method'] = 'SPIDER'
    spider_corr_sub = spider_corr_sub.sort_values('correlation', ascending=False)
    corr_sub = pd.concat([spider_corr_sub, sdm_corr_sub])
    # corr_sub = corr_sub.query('correlation > 0.3')
    import seaborn as sns
    plt.figure(figsize=(5, 10))
    sns.barplot(data=corr_sub, y='celltype', x='correlation', hue='method', hue_order=['SPIDER', 'SpatialDM'], palette={'SpatialDM':'#313695', 'SPIDER': '#D53E4F'}, orient='h')
    plt.xticks(rotation=90)
    
def plot_sdm_per_ct(adata, sdm_adata, idata):
    import matplotlib.pyplot as plt
    st_decon = adata.obsm['deconvolution']
    st_decon_org = st_decon.copy()
    decon_corr = pd.concat([st_decon, sdm_adata.to_df()], axis=1).corr()[sdm_adata.var_names].loc[st_decon.columns].T
    decon_corr_sup = decon_corr.loc[decon_corr.idxmax()] 
    decon_corr = decon_corr[decon_corr_sup.columns[np.argsort(-decon_corr_sup.to_numpy().diagonal())]]
    decon_corr_sup = decon_corr.loc[decon_corr.idxmax()] 
    plt.rcParams.update({'font.size': 25})
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

    marker_score = sdm_adata.to_df()[decon_corr_sup.index.unique()]
    marker_score = pd.concat([sdm_adata.obs[['x', 'y']], marker_score], axis=1)
    for i in range(len(decon_corr_sup.index)):
        plt.subplot(2, len(celtypes), base)
        im=plt.scatter(marker_score['x'],marker_score['y'], c=marker_score[decon_corr_sup.index[i]].to_numpy(), s=50, cmap='plasma')
        plt.colorbar(im,fraction=0.046, pad=0.04)
        plt.axis('equal')
        plt.axis('off')
        plt.title(f'{decon_corr_sup.index[i]} \n corr={"%.3f" % corrs[i]}')
        base += 1

    plt.tight_layout()
    plt.show()
    
def plot_spider_per_ct(idata, adata):
    import anndata
    import matplotlib.pyplot as plt
    idata = idata[:, idata.var['is_svi'] == 1]
    op.util.interaction_spot_interface(idata)
    adata = adata[adata.obs_names.isin(idata.uns['cell_score'].index)]
    adata_lri = anndata.AnnData(idata.uns['cell_score'].loc[adata.obs_names])
    adata_lri.obs=adata.obs
    adata_lri.obsm['spatial'] = adata.obsm['spatial']
    adata_lri = adata_lri[:, adata_lri.var_names[idata.var['is_svi'].to_numpy().astype(bool)]]
    st_decon = adata.obsm['deconvolution']
    st_decon_org = st_decon.copy()
    decon_corr = pd.concat([st_decon, adata_lri.to_df()], axis=1).corr()[adata_lri.var_names].loc[st_decon.columns].T
    decon_corr_sup = decon_corr.loc[decon_corr.idxmax()] 
    decon_corr = decon_corr[decon_corr_sup.columns[np.argsort(-decon_corr_sup.to_numpy().diagonal())]]
    decon_corr_sup = decon_corr.loc[decon_corr.idxmax()] 
    plt.rcParams.update({'font.size': 25})
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

    marker_score = adata_lri.to_df()[decon_corr_sup.index.unique()]
    marker_score = pd.concat([adata_lri.obs[['x', 'y']], marker_score], axis=1)
    for i in range(len(decon_corr_sup.index)):
        plt.subplot(2, len(celtypes), base)
        im=plt.scatter(marker_score['x'],marker_score['y'], c=marker_score[decon_corr_sup.index[i]].to_numpy(), s=50, cmap='plasma')
        plt.colorbar(im,fraction=0.046, pad=0.04)
        plt.axis('equal')
        plt.axis('off')
        plt.title(f'{decon_corr_sup.index[i]} \n corr={"%.3f" % corrs[i]}')
        base += 1

    plt.tight_layout()
    plt.show()
    
def plot_sdm_pattern(adata, sdm_patterns):
    plt.figure(figsize=(15,2))
    import spatialdm.plottings as pl
    plt.rcParams.update({'font.size': 13})
    for i in range(5):
        plt.subplot(1, 5, i + 1)
        plt.scatter(adata.obsm['spatial'][:,0], adata.obsm['spatial'][:,1], c=sdm_patterns[str(i)], s=25, cmap='plasma',);
        plt.axis('equal')
        plt.axis('off')
        pl.plt_util('Pattern {}'.format(i))
    plt.show()
    
def ct_pattern_corr(idata, adata, sdm_patterns):
    adata, adata_lri, adata_pattern = op.cell_transform(idata, adata, label=adata.uns['cluster_key'])
    ldf = pd.get_dummies(adata.obs['cluster'])
    plt.figure(figsize=(3, 4))
    plt.rcParams.update({'font.size': 14})
    plt.subplot(2, 1, 1)
    sns.heatmap(np.corrcoef(ldf.to_numpy(), adata.obsm['interaction_pattern'],  rowvar=False)[:4, 4:], vmin=-1, vmax=1, cmap='coolwarm', yticklabels=ldf.columns.str.replace('_', ' '),)
    plt.ylabel('SPIDER\nSVI pattern')
    plt.title('Correlations')
    plt.subplot(2, 1, 2)
    ldf = pd.get_dummies(adata.obs['cluster'])
    sns.heatmap(np.corrcoef(ldf.to_numpy(), sdm_patterns, rowvar=False)[:4, 4:], vmin=-1, vmax=1, cmap='coolwarm', yticklabels=ldf.columns.str.replace('_', ' '),)
    plt.ylabel('SpatialDM\nLRI pattern')
    plt.show()
    
def ct_cross_corr(idata):
    import networkx as nx
    import anndata
    from matplotlib import pyplot as plt
    pattern_idata = anndata.AnnData(idata.obsm['pattern_score'])
    pattern_idata.obs_names = idata.obs_names
    pattern_idata.obs = idata.obs
    pattern_idata.obs['A_label'] = pattern_idata.obs['A_label'].str.replace('_', ' ')
    pattern_idata.obs['B_label'] = pattern_idata.obs['B_label'].str.replace('_', ' ')
    label_1 = pattern_idata.obs['A_label'].astype(str) + '_' + pattern_idata.obs['B_label'].astype(str).to_numpy()
    label_2 = pattern_idata.obs['B_label'].astype(str) + '_' + pattern_idata.obs['A_label'].astype(str).to_numpy()
    pick = pattern_idata.obs[['label_1', 'label_2']].astype(int).idxmax(axis=1).to_numpy()
    text_label = [label_1[i] if x=='label_1' else label_2[i] for i,x in enumerate(pick)]
    pattern_idata.obs['label'] = text_label
    pattern_idata.obs
    pattern_idata.var_names = ['Pattern {}'.format(x) for x in range(pattern_idata.shape[1])]
    plt.figure(figsize=(len(pattern_idata.var_names)*2, 1.5))
    count = 1
    for pattern_id in pattern_idata.var_names:
        df = pd.concat([pattern_idata.to_df()[pattern_id],  pattern_idata.obs], axis=1)
        df[pattern_id] = (df[pattern_id] - df[pattern_id].min()) / (df[pattern_id].max() - df[pattern_id].min())
        input_df = df.groupby('label')[pattern_id].mean().reset_index(drop=False)
        input_df[['source', 'target']] = input_df['label'].str.split('_', expand=True).to_numpy() 
        links_df = input_df[['source', 'target', pattern_id]]
        links_df.columns = ['source', 'target', 'value']

        # reshape to a adjacency matrix/
        links_df = links_df.pivot(index='source', columns='target', values='value').replace(np.nan, 0)
        import seaborn as sns
        import matplotlib.pyplot as plt
        mask = np.triu(np.ones(links_df.shape),1).astype(bool)
        plt.subplot(1, len(pattern_idata.var_names), count)
        sns.heatmap(links_df, cmap='plasma', cbar=True, mask=mask, fmt='.2f', linewidths=0.5, annot_kws={'size': 15})
        # remove y ticks except the first subplot
        if count != 1:
            plt.yticks([])
        plt.ylabel('')
        plt.xlabel('')
        plt.title(pattern_id)
        count += 1
    plt.tight_layout() 
    plt.show()  
    
def get_interaction_type(svi_df, is_human):
    if is_human:
        cellchat = pd.read_csv('../CellChatDB/1_human/interaction_input_CellChatDB.csv', index_col=0)
        cellchat.columns = [x.replace('interaction.','') for x in cellchat.columns]
        cellchat = cellchat[['ligand', 'receptor', 'annotation']].reset_index()
        cellchat[['Ligand0', 'Receptor0', 'Receptor1', 'Receptor2', 'Receptor3']] = cellchat['index'].str.split('_', expand=True)
        cellchat.set_index('index', inplace=True)
        expand_cellchat = []
        addition = []
        for x in cellchat.to_numpy():
            _, r, anno, l, r1, r2, r3, r4 = x
            if r4:
                addition.append([l, r4, anno])
            if r3:
                addition.append([l, r3, anno])
            if r2:
                addition.append([l, r2, anno])
            expand_cellchat.append([l, r1, anno])
        expand_cellchat = pd.DataFrame(expand_cellchat, columns=['ligand', 'receptor', 'annotation'])
        expand_cellchat.index = expand_cellchat['ligand'] + '_' + expand_cellchat['receptor']
    else:
        cellchat = pd.read_csv('../CellChatDB/2_mouse/interaction_input_CellChatDB.csv', index_col=0)
        cellchat = cellchat[['ligand', 'receptor', 'annotation']].reset_index()
        cellchat[['Ligand0', 'Receptor0', 'Receptor1', 'Receptor2', 'Receptor3']] = cellchat['index'].str.split('_', expand=True)
        cellchat.set_index('index', inplace=True)
        expand_cellchat = []
        addition = []
        for x in cellchat.to_numpy():
            _, r, anno, l, r1, r2, r3, r4 = x
            if r4:
                addition.append([l, r4, anno])
            if r3:
                addition.append([l, r3, anno])
            if r2:
                addition.append([l, r2, anno])
            expand_cellchat.append([l, r1, anno])
        expand_cellchat = pd.DataFrame(expand_cellchat, columns=['ligand', 'receptor', 'annotation'])
        expand_cellchat['ligand'] = expand_cellchat['ligand'].str.capitalize()
        expand_cellchat['receptor'] = expand_cellchat['receptor'].str.capitalize()
        expand_cellchat.index = expand_cellchat['ligand'] + '_' + expand_cellchat['receptor']
        
    joint = svi_df.join(expand_cellchat[['annotation']], how='left')
    return joint, joint['annotation'].value_counts()

def plot_interaction_type(spider_counts_df, spatialdm_counts_df):
    counts_df = pd.concat([spider_counts_df, spatialdm_counts_df], axis=1)
    counts_df.columns = ['SPIDER', 'SpatialDM']
    counts_df = counts_df.div(counts_df.sum(axis=0), axis=1)
    melted_counts = counts_df.reset_index().melt(value_name='Percentage', var_name='Method', id_vars='index')
    plt.figure(figsize=(10, 5))
    sns.barplot(data=melted_counts, palette={'SpatialDM':'#313695', 'SPIDER': '#D53E4F'}, x='index', y='Percentage', hue='Method')

def plot_interaciton_hetro_homo(idata):
    plt.figure(figsize=(8, 5))
    svi_count_df = idata.to_df()[idata.var.query('is_svi==1').index]
    # normalize each column
    svi_count_df = (svi_count_df - svi_count_df.min()) / (svi_count_df.max() - svi_count_df.min())
    svi_count_df['label'] = (idata.obs['A_label'] == idata.obs['B_label']).astype(str).replace('True', 'Homotypic').replace('False', 'Heterotypic')
    sig_type = svi_count_df.groupby('label').median()
    sig_frac = (sig_type / sig_type.sum()).T
    sns.boxplot(data=sig_frac.melt(), x='label', y='value', palette={'Heterotypic':'#F7994A', 'Homotypic': '#7F46C0'})
    plt.ylabel('Fraction of total signaling')
    plt.xlabel('')
    
def plot_lris_homo(idata, sig_frac_homo):
    plt.figure(figsize=(15,5))
    plt.rcParams.update({'font.size': 13})
    fig_count = 1
    for (i, value) in sig_frac_homo[:5].reset_index()[['index', 'Homotypic']].to_numpy():
        plt.subplot(2, 5, fig_count)
        fig_count += 1
        plt.scatter(idata.obsm['spatial'][:,0], idata.obsm['spatial'][:,1], c=idata.to_df()[i], s=10, cmap='plasma',);
        plt.axis('equal')
        plt.axis('off')
        plt.title(f'{i}\nfrac={"%.3f" % value}', y=0.95)
    plt.show()
    
def plot_lris_hete(idata, sig_frac_heter):
    plt.figure(figsize=(15,5))
    plt.rcParams.update({'font.size': 13})
    fig_count = 1
    for (i, value) in sig_frac_heter[:5].reset_index()[['index', 'Heterotypic']].to_numpy():
        plt.subplot(2, 5, fig_count)
        fig_count += 1
        plt.scatter(idata.obsm['spatial'][:,0], idata.obsm['spatial'][:,1], c=idata.to_df()[i], s=10, cmap='plasma',);
        plt.axis('equal')
        plt.axis('off')
        plt.title(f'{i}\nfrac={"%.3f" % value}', y=0.95)
    plt.show()
    
    
def plot_percentage_hete_homo(idata):    
    svi_count_df = idata.to_df()[idata.var.query('is_svi==1').index]
    svi_count_df = (svi_count_df - svi_count_df.min()) / (svi_count_df.max() - svi_count_df.min())
    svi_count_df['label'] = (idata.obs['A_label'] == idata.obs['B_label']).astype(str).replace('True', 'Homotypic').replace('False', 'Heterotypic')
    sig_type = svi_count_df.groupby('label').median()
    sig_frac = (sig_type / sig_type.sum()).T
    sig_frac['type'] = 'Homotypic'
    sig_frac.loc[sig_frac['Heterotypic'] > sig_frac['Homotypic'], 'type'] = 'Heterotypic'
    
    sig_frac['diff'] = sig_frac['Heterotypic'] - sig_frac['Homotypic']
    sig_frac_heter = sig_frac.sort_values('diff', ascending=False)[:10]
    sig_frac_homo = sig_frac.sort_values('diff', ascending=True)[:10]
    sig_frac_selected = pd.concat([sig_frac_homo, sig_frac_heter])
    sig_frac_selected[['Heterotypic', 'Homotypic']].plot(kind='bar', stacked=True, color=['#F7994A', '#7F46C0'], figsize=(9, 5))
    return sig_frac, sig_frac_homo, sig_frac_heter


def ct_chord(idata, svi):
    import holoviews as hv
    from holoviews import opts, dim
    df = pd.concat([idata.to_df()[svi],  idata.obs], axis=1)
    df['direction'] = np.array(idata[:, idata.var_names==svi].layers['direction']).flatten()
    df['arrow'] = df.apply(lambda x: f'{x.A_label}->{x.B_label}' if x.direction==0 else f'{x.B_label}->{x.A_label}', axis=1)
    input_df = df.groupby('arrow')[svi].median().reset_index(drop=False).query(f'{svi} > 0')
    input_df[['source', 'target']] = input_df['arrow'].str.split('->', expand=True).to_numpy() 
    node_df = pd.DataFrame(index=np.unique(np.concatenate((input_df['source'].to_numpy(), input_df['target'].to_numpy()))))
    node_df['group'] = 1
    links_df = input_df[['source', 'target', svi]]
    links_df.columns = ['source', 'target', 'value']

    hv.extension('bokeh')
    hv.output(size=200)

    nodes = hv.Dataset(node_df, 'index')
    def rotate_label(plot, element):
        angles = plot.handles['text_1_source'].data['angle']
        characters = np.array(plot.handles['text_1_source'].data['text'])
        plot.handles['text_1_source'].data['text'] = np.array([x + ' ' * int(len(x)) if x in characters[np.where((angles < -1.5707963267949) | (angles > 1.5707963267949))] else x for x in plot.handles['text_1_source'].data['text']])
        new_text = []
        for x in plot.handles['text_1_source'].data['text']:
            new = x
            if x in characters[np.where((angles > -1.5707963267949) & (angles < 1.5707963267949))]:
                new = ' ' * int(len(x)) + new
            if x in characters[np.where((angles > -1.5707963267949) & (angles < 0))]:
                # new = '\n' + new
                new = '\n\n' + new
            else:
                # new = new + '\n'
                new = new + '\n\n'
            new_text.append(new)
        plot.handles['text_1_source'].data['text'] = np.array(new_text)
        angles[np.where((angles > 1.5707963267949))] = 0 # left top
        angles[np.where((angles < -1.5707963267949))] = 0  # left bottom
        angles[np.where((angles < 1.5707963267949)  & (angles > 0) )] = 0 # right top
        angles[np.where((angles > -1.5707963267949) & (angles < 0))] = 0 # right bottom
        # angles = np.zeros(len(angles))
        plot.handles['text_1_glyph'].text_align = "center"
    chord = hv.Chord((links_df, nodes))
    chord.opts(
        opts.Chord(cmap='Set2', edge_cmap='Set2', edge_color=dim('source').str(), 
                labels='index', node_color=dim('index').str(),   hooks=[rotate_label], label_text_font_size='17pt', title=f'{svi}'),
    )
    return chord

def mouse_brain_pattern_corr(idata, adata):
    idata = idata[:, idata.var['is_svi']==1]
    adata, adata_lri, adata_pattern = op.cell_transform(idata, adata, label=adata.uns['cluster_key'])
    adata_lri.obsm['spatial'] = adata[adata_lri.obs_names].obsm['spatial']
    adata_pattern.obsm['spatial'] = adata[adata_lri.obs_names].obsm['spatial']

    region = adata.obsm['deconvolution_results'].loc[adata_pattern.obs_names]
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
        'Ependymal': 'Ependymal',
        'CA1_CA2_CA3_Subiculum': "CA1/CA2/CA3 Subiculum",
        'Subiculum_Entorhinal_cl2': "Subiculum Entorhinal cl2"
    }
    for x in useful_df:
        patterns.append(int(x.idxmax()))
        corrs.append(x.max())
        cluster_id.append(x.name)
        cluster.append(cluster_id_map[x.name])
    patterns = np.array(patterns)[np.argsort(corrs)[::-1]]
    cluster_id = np.array(cluster_id)[np.argsort(corrs)[::-1]]
    cluster = np.array(cluster)[np.argsort(corrs)[::-1]]
    corrs = np.array(corrs)[np.argsort(corrs)[::-1]]
    
    useful_df = []
    adata_lri.obsm['spatial'] = adata[adata_lri.obs_names].obsm['spatial']
    adata_pattern.obsm['spatial'] = adata[adata_lri.obs_names].obsm['spatial']

    plt.figure(figsize=(7*4, 5*len(patterns)))
    plt.rcParams['font.size'] = 25
    base = 1
    for i in range(len(patterns)):
        corr =  corrs[i]
        cluster_name = cluster[i]
        cluster_i = cluster_id[i]

        plt.subplot(len(patterns), 6, base)

        # plt.scatter(adata.obsm['spatial'][:,0],adata.obsm['spatial'][:,1], c=adata.obsm['deconvolution_results'][cluster_i], s=3, cmap='plasma', edgecolors='none')
        to_plot = pd.DataFrame(adata.obsm['spatial'], columns = ['x', 'y'])
        to_plot['cluster'] = adata.obsm['deconvolution_results'][cluster_i].to_numpy()
        to_plot = to_plot.sort_values('cluster', ascending=True)
        plt.scatter(to_plot['x'].to_numpy(),to_plot['y'].to_numpy(), c=to_plot['cluster'].to_numpy(), s=5, cmap='plasma', vmin=0, edgecolors='none')

        # plt.scatter(adata[adata.obs['cluster']==cluster_id].obsm['spatial'][:,0],adata[adata.obs['cluster']==cluster_id].obsm['spatial'][:,1], c=adata[adata.obs['cluster']==cluster_name].obs['cluster']==cluster_name, s=5, cmap='plasma', vmin=0, edgecolors='none')
        plt.axis('equal')
        plt.axis('off')
        plt.title(cluster_name.replace(" ", "\n", 1), y=0.92)
        base += 1
        
        plt.subplot(len(patterns), 6, base)
        im=plt.scatter(adata_pattern.obsm['spatial'][:,0], adata_pattern.obsm['spatial'][:,1], c=adata_pattern.to_df()[str(patterns[i])], s=3, cmap='plasma', edgecolors='none')
        plt.colorbar(im,fraction=0.046, pad=0.04)
        plt.axis('equal')
        plt.axis('off')
        plt.title(f'Pattern {patterns[i]} - {np.sum(idata.var.label == patterns[i])} SVIs\ncorr={"%.3f" % corr}', y=0.9)
        base += 1
        
        svis = idata.var[idata.var.label == patterns[i]].sort_values(f'pattern_correlation_{patterns[i]}', ascending=False)
        svi_to_plot = svis.index.to_numpy()[:4]
        print(svi_to_plot)
        memberships = svis[f'pattern_correlation_{patterns[i]}'].to_numpy()[:4]
        for j in range(4):
            if j < len(svi_to_plot):
                plt.subplot(len(patterns), 6, base)
                im=plt.scatter(adata_lri.obsm['spatial'][:,0], adata_lri.obsm['spatial'][:,1], c=adata_lri.to_df()[svi_to_plot[j]], s=3, cmap='plasma', edgecolors='none',linewidths=0)
                plt.axis('equal')
                plt.axis('off')
                plt.colorbar(im,fraction=0.046, pad=0.04)
                plt.title(f'{svi_to_plot[j]}\n corr={"%.3f" % memberships[j]}', y=0.9)
            base += 1
            
def saw_celltype_corr(idata, adata):
    import matplotlib.pyplot as plt 
    adata, adata_lri, adata_pattern = op.cell_transform(idata[:, idata.var.query('is_svi==1').index], adata, label=adata.uns['cluster_key'])
    st_decon = pd.get_dummies(adata_lri.obs['cell_type'])
    st_decon_org = st_decon.copy()
    decon_corr = pd.concat([st_decon, adata_lri.to_df()], axis=1).corr()[adata_lri.var_names].loc[st_decon.columns].T
    decon_corr_sup = decon_corr.loc[decon_corr.idxmax()] 
    decon_corr = decon_corr[decon_corr_sup.columns[np.argsort(-decon_corr_sup.to_numpy().diagonal())]]
    decon_corr_sup = decon_corr.loc[decon_corr.idxmax()] 
    plt.rcParams.update({'font.size': 23})
    st_decon = st_decon_org.copy()
    celtypes = decon_corr_sup.columns
    st_decon = pd.concat([idata.uns['cell_meta'][['x', 'y']], st_decon], axis=1)
    plt.figure(figsize=(len(celtypes)*4, 8))
    base = 1
    for i in celtypes:
        plt.subplot(2, len(celtypes), base)
        # sns.scatterplot(st_decon, x='x', y='y', hue=i, s=10, linewidth=0)
        im=plt.scatter(st_decon['x'],st_decon['y'], c=st_decon[i], s=10, cmap='plasma')
        plt.colorbar(im,fraction=0.046, pad=0.04)

        plt.axis('equal')
        plt.axis('off')
        plt.title(i.replace(' ', '\n'))
        base += 1
    corrs = decon_corr_sup.to_numpy().diagonal()
    marker_score = adata_lri.to_df()[decon_corr_sup.index.unique()]
    marker_score = pd.concat([adata_lri.obs[['x', 'y']], marker_score], axis=1)
    for i in range(len(decon_corr_sup.index)):
        plt.subplot(2, len(celtypes), base)
        im=plt.scatter(marker_score['x'],marker_score['y'], c=marker_score[decon_corr_sup.index[i]].to_numpy(), s=10, cmap='plasma')
        plt.colorbar(im,fraction=0.046, pad=0.04)
        plt.axis('equal')
        plt.axis('off')
        plt.title(f'{decon_corr_sup.index[i]}\ncorr={"%.3f" % corrs[i]}')
        base += 1
    plt.tight_layout()

def multiple_leiden(adata_lri, adata):
    import seaborn as sns
    import matplotlib.pyplot as plt 
    plt.figure(figsize=(7, 2))
    plt.subplot(1, 4, 1)
    import matplotlib.pyplot as plt
    to_plot = pd.DataFrame(adata.obs['cell_type'])
    to_plot[['x','y']] = adata.obsm['spatial']
    import seaborn as sns
    custom_palette = sns.color_palette("husl",to_plot.cell_type.nunique())
    sns.scatterplot(data=to_plot, x='x', y='y', hue='cell_type', palette=custom_palette, s=3, linewidth=0, legend=False)
    plt.axis('equal')
    plt.axis('off')

    plt.subplot(1, 4, 2)
    to_plot = pd.DataFrame(adata_lri.obs['leiden_res1'])
    to_plot[['x','y']] = adata_lri.obsm['spatial']
    custom_palette = sns.color_palette("husl",to_plot['leiden_res1'].nunique())
    sns.scatterplot(data=to_plot, x='x', y='y', hue='leiden_res1', palette=custom_palette, s=3, linewidth=0, legend=True)
    plt.legend(bbox_to_anchor=(1, 0), loc=1, borderaxespad=0., ncols=5, frameon=False, fontsize=5)
    plt.title('Leiden\nresolution=1')
    plt.axis('equal')
    plt.axis('off')

    plt.subplot(1, 4, 3)
    to_plot = pd.DataFrame(adata_lri.obs['leiden_res0.5'])
    to_plot[['x','y']] = adata_lri.obsm['spatial']
    custom_palette = sns.color_palette("husl",to_plot['leiden_res0.5'].nunique())
    sns.scatterplot(data=to_plot, x='x', y='y', hue='leiden_res0.5', palette=custom_palette, s=3, linewidth=0, legend=True)
    plt.title('Leiden\nresolution=0.5')
    plt.legend(bbox_to_anchor=(1, 0), loc=1, borderaxespad=0., ncols=5, frameon=False, fontsize=5)
    plt.axis('equal')
    plt.axis('off')

    plt.subplot(1, 4, 4)
    to_plot = pd.DataFrame(adata_lri.obs['leiden_res0.1'])
    to_plot[['x','y']] = adata_lri.obsm['spatial']
    custom_palette = sns.color_palette("husl",to_plot['leiden_res0.1'].nunique())
    sns.scatterplot(data=to_plot, x='x', y='y', hue='leiden_res0.1', palette=custom_palette, s=3, linewidth=0, legend=True)
    plt.legend(bbox_to_anchor=(1, 0), loc=1, borderaxespad=0., ncols=5, frameon=False, fontsize=5)
    plt.title('Leiden\nresolution=0.1')
    plt.axis('equal')
    plt.axis('off')
    
def multiple_leiden_oligo(adata_lri, adata):
    import seaborn as sns
    import matplotlib.pyplot as plt 
    plt.figure(figsize=(7, 2))
    plt.subplot(1, 4, 1)
    import matplotlib.pyplot as plt
    to_plot = pd.DataFrame(adata.obs['cell_type'])
    to_plot[['x','y']] = adata.obsm['spatial']
    import seaborn as sns
    custom_palette = sns.color_palette("husl",to_plot.cell_type.nunique())
    sns.scatterplot(data=to_plot, x='x', y='y', hue='cell_type', palette=custom_palette, s=3, linewidth=0, legend=False)
    # move legend to right
    # plt.legend(bbox_to_anchor=(0, 0), loc=1, borderaxespad=0., ncols=1, frameon=False, fontsize=3)
    plt.axis('equal')
    plt.axis('off')

    plt.subplot(1, 4, 2)
    to_plot = pd.DataFrame(adata_lri.obs['leiden_res1'])
    to_plot[['x','y']] = adata_lri.obsm['spatial']
    to_plot['leiden_res1_oligo'] = to_plot['leiden_res1']
    to_plot.loc[adata_lri.obs[adata_lri.obs['cell_type'] != 'Oligodendrocyte'].index, 'leiden_res1_oligo'] = np.nan
    sns.scatterplot(data=to_plot, x='x', y='y', s=15, linewidth=0, legend=False, color='grey')
    custom_palette = sns.color_palette("husl",to_plot['leiden_res1'].nunique())
    sns.scatterplot(data=to_plot, x='x', y='y', hue='leiden_res1_oligo', s=5, linewidth=0, legend=True, palette=custom_palette)
    plt.title('Oligo\nresolution=1', fontsize=8)
    plt.legend(bbox_to_anchor=(1, 0), loc=1, borderaxespad=0., ncols=5, frameon=False, fontsize=5)
    plt.axis('equal')
    plt.axis('off')

    plt.subplot(1, 4, 3)
    to_plot = pd.DataFrame(adata_lri.obs['leiden_res0.5'])
    to_plot[['x','y']] = adata_lri.obsm['spatial']
    to_plot['leiden_res0.5_oligo'] = to_plot['leiden_res0.5']
    to_plot.loc[adata_lri.obs[adata_lri.obs['cell_type'] != 'Oligodendrocyte'].index, 'leiden_res0.5_oligo'] = np.nan
    sns.scatterplot(data=to_plot, x='x', y='y', s=15, linewidth=0, legend=False, color='grey')
    custom_palette = sns.color_palette("husl",to_plot['leiden_res0.5'].nunique())
    sns.scatterplot(data=to_plot, x='x', y='y', hue='leiden_res0.5_oligo', s=5, linewidth=0, legend=True, palette=custom_palette)
    plt.title('Oligo\nresolution=0.5', fontsize=8)
    plt.legend(bbox_to_anchor=(1, 0), loc=1, borderaxespad=0., ncols=5, frameon=False, fontsize=5)
    plt.axis('equal')
    plt.axis('off')

    plt.subplot(1, 4, 4)
    to_plot = pd.DataFrame(adata_lri.obs['leiden_res0.1'])
    to_plot[['x','y']] = adata_lri.obsm['spatial']
    to_plot['leiden_res0.1_oligo'] = to_plot['leiden_res0.1']
    to_plot.loc[adata_lri.obs[adata_lri.obs['cell_type'] != 'Oligodendrocyte'].index, 'leiden_res0.1_oligo'] = np.nan
    sns.scatterplot(data=to_plot, x='x', y='y', s=15, linewidth=0, legend=False, color='grey')
    custom_palette = sns.color_palette("husl",to_plot['leiden_res0.1'].nunique())
    sns.scatterplot(data=to_plot, x='x', y='y', hue='leiden_res0.1_oligo', s=5, linewidth=0, legend=True, palette=custom_palette)
    plt.legend(bbox_to_anchor=(1, 0), loc=1, borderaxespad=0., ncols=5, frameon=False, fontsize=5)
    plt.title('Oligo\nresolution=0.1', fontsize=8)
    plt.axis('equal')
    plt.axis('off')
    
def oligo_marker(adata_lri):
    import scanpy as sc
    adata_lri.obs['leiden_res0.1_oligo'] = adata_lri.obs['leiden_res0.1']
    adata_lri.obs.loc[adata_lri.obs[adata_lri.obs['cell_type'] != 'Oligodendrocyte'].index, 'leiden_res0.1_oligo'] = np.nan
    adata_lri_valuecounts = adata_lri.obs['leiden_res0.1_oligo'].value_counts()
    labels = adata_lri_valuecounts.index[adata_lri_valuecounts > 1].to_numpy()
    adata_lri_subset = adata_lri[adata_lri.obs['leiden_res0.1_oligo'].isin(labels)]
    sc.tl.rank_genes_groups(adata_lri_subset, groupby='leiden_res0.1_oligo')
    sc.set_figure_params(dpi=300, fontsize=10, figsize=(10, 10))
    sc.pl.rank_genes_groups_dotplot(adata_lri_subset, standard_scale='var', n_genes=1, show=False, dendrogram=False, color_map='plasma', swap_axes=True, title='Oligo\nresolution=0.1')

    adata_lri.obs['leiden_res0.5_oligo'] = adata_lri.obs['leiden_res0.5']
    adata_lri.obs.loc[adata_lri.obs[adata_lri.obs['cell_type'] != 'Oligodendrocyte'].index, 'leiden_res0.5_oligo'] = np.nan
    adata_lri_valuecounts = adata_lri.obs['leiden_res0.5_oligo'].value_counts()
    labels = adata_lri_valuecounts.index[adata_lri_valuecounts > 1].to_numpy()
    adata_lri_subset = adata_lri[adata_lri.obs['leiden_res0.5_oligo'].isin(labels)]
    sc.tl.rank_genes_groups(adata_lri_subset, groupby='leiden_res0.5_oligo')
    sc.set_figure_params(dpi=300, fontsize=10, figsize=(10, 10))
    sc.pl.rank_genes_groups_dotplot(adata_lri_subset, standard_scale='var', n_genes=1, show=False, dendrogram=False, color_map='plasma', swap_axes=True, title='Oligo\nresolution=0.5')


    adata_lri.obs['leiden_res1_oligo'] = adata_lri.obs['leiden_res1'].to_numpy()
    adata_lri.obs.loc[adata_lri.obs[adata_lri.obs['cell_type'] != 'Oligodendrocyte'].index, 'leiden_res1_oligo'] = np.nan
    adata_lri_valuecounts = adata_lri.obs['leiden_res1_oligo'].value_counts()
    labels = adata_lri_valuecounts.index[adata_lri_valuecounts > 1].to_numpy()
    adata_lri_subset = adata_lri[adata_lri.obs['leiden_res1_oligo'].isin(labels)]
    sc.tl.rank_genes_groups(adata_lri_subset, groupby='leiden_res1_oligo')
    sc.set_figure_params(dpi=300, fontsize=10, figsize=(10, 10))
    sc.pl.rank_genes_groups_dotplot(adata_lri_subset, standard_scale='var', n_genes=1, show=False, dendrogram=False, color_map='plasma', swap_axes=True, title='Oligo\nresolution=1')
    
    
def mouse_brain_ct_comp(adata_pattern_subset):
    dfs = []
    for i in adata_pattern_subset.obs['leiden'].unique():
        dfs.append(adata_pattern_subset.obs['cluster'].loc[(adata_pattern_subset.obs[adata_pattern_subset.obs['leiden']==i]).index].value_counts())
    dfs_merged = pd.concat(dfs, axis=1)
    dfs_merged.columns=adata_pattern_subset.obs['leiden'].unique()
    dfs_merged = (dfs_merged / dfs_merged.sum()).T.reset_index().rename(columns={'index':'interface cluster'}).sort_values('Astrocytes')
    # get largest and second largest in each row of df
    subset=dfs_merged[(dfs_merged.set_index('interface cluster').apply(lambda row: row.nlargest(2).values[-1],axis=1) >= 0.2).to_numpy()]
    dfs_merged = dfs_merged[(dfs_merged.set_index('interface cluster').max(axis=1)>=0.3).to_numpy()]
    custom_palette = sns.color_palette("husl",adata_pattern_subset.obs['cluster'].nunique())
    plt.rcParams['font.size'] = 12
    plt.figure(figsize=(20, 15))
    # g = dfs_merged.plot.barh(x='interface cluster',
    g = dfs_merged.set_index('interface cluster').reset_index().plot.barh(x='interface cluster',
            stacked=True, legend=False, color=custom_palette, figsize=(5,2), width=0.9)
    # g = dfs_merged.set_index('interface cluster').loc[np.array([7, 1, 8, 9, 13, 16, 11, 3, 14, 19]).astype(str)].reset_index().plot.barh(x='interface cluster',
    #         stacked=True, legend=False, color=custom_palette, figsize=(6,3), width=0.9)
    g.spines[['right', 'top']].set_visible(False)
    plt.ylabel('Clusters from\nSVI patterns')
    plt.title('Celltype compositions')
    plt.legend(bbox_to_anchor=(0.98, -0.1), ncols=3, fontsize=9, frameon=False)
    labels = []
    for j in dfs_merged.columns:
        if j != 'interface cluster':
            for i in dfs_merged.index:
                label = j
                labels.append(label)
    for label, p in zip(labels, g.patches):
        width, height = p.get_width(), p.get_height()
        x, y = p.get_xy() 
        if width>=0.3:
            if label == 'Subiculum_Entorhinal_cl2':
                g.text(x+width/2-0.1,  y+height/2-0.01, 
                label, # '{:.0f} %'.format(height), 
                horizontalalignment='center', 
                verticalalignment='center',
                # rotation=90,
                )
            else:
                g.text(x+width/2-0.02,  y+height/2-0.01, 
                label, # '{:.0f} %'.format(height), 
                horizontalalignment='center', 
                verticalalignment='center')
    