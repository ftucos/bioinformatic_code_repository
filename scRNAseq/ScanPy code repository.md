# ScanPy code repository

## Extract complete PCs contributions by different genes

```python
sc.pl.pca_overview(adata)

# or extract as a full list
components = pd.DataFrame(adata.varm['PCs'][:,0:2])
components.index = adata.var.index
components.columns = ['PC1', 'PC2']
components.sort_values('PC1')

```

## Select highly variable genes with triku in place of variance/mean relation: 

```python
# Select features with triku

import triku as tk

sc.pp.pca(adata)
sc.pp.neighbors(adata)

tk.tl.triku(adata, n_features=5000)
```

### Add PCA components to .obs in order to be explored and plotted to debug low dimensionality embedding

```python
for i in range(30):
	adata.obs['PCA_' + str(i)] = adata.obsm['X_pca'][:,i]
  
  sc.pl.scatter(adata, basis='umap', color='PCA_20')

```

or 





### Plot venn diagram intersection

```python
from matplotlib_venn import venn3

venn3([set(wc),set(tt),set(tt_ov)], ('Wilcox','T-test','T-test_ov') )
plt.show()
```

### Plot Heatmap with standard deviation scale

```python
sc.settings.figdir = "output/Markers/NUMB_WT/"
sc.tl.rank_genes_groups(adata_PID_wt, groupby='Clusters_simplified', key_added='rank_genes_clusters', method='wilcoxon')
sc.pl.rank_genes_groups_heatmap(sc.pp.scale(adata_PID_wt, max_value=10, copy=True),
								key='rank_genes_clusters', show_gene_labels=True,
								use_raw=False, n_genes = 10,
								cmap='RdBu_r', vmin=-3, vmax=3, save="_top_10-NUMB_WT-PID_pathways_by_cluster.png") # _r reverses the colormap direction
```

