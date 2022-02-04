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

