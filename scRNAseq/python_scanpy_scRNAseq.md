# scRNAseq analysis in python

Install all the  required packages. Whenever possible use conda in place of pip since it will deal better with dependencies versions compatibilities.

```shell	
conda install -c bioconda conda-forge numpy pandas matplotlib seaborn scikit-learn  statsmodels numba pytables python-igraph leidenalg louvain scipy anndata scanpy dca scikit-learn-intelex hnswlib pybind11 notebook jupyterlab
conda install keras=2.4
pip install git+https://github.com/theislab/scvelo
```

## Import raw matrices into an andata object

**Define import function**

Import function to load data from star output directory

```python
def buildAnndataFromStar(path):
    """Generate an anndata object from the STAR aligner output folder"""
    path=path
    # Load Read Counts
    X = sc.read_mtx(path+'/Gene/raw/matrix.mtx')
    
    # Transpose counts matrix to have Cells as rows and Genes as cols as expected by AnnData objects
    X = X.X.transpose()
    
    # This matrix is organized as a sparse matrix with Row, Colm and 3 values colums for 
    # Spliced, Unspliced and ambigous reads
    mtx = np.loadtxt(path+'/Velocyto/raw/matrix.mtx', skiprows=3, delimiter=' ')
    # Extract sparse matrix shape informations from the third row
    shape = np.loadtxt(path+'/Velocyto/raw/matrix.mtx', skiprows=2, max_rows = 1 ,delimiter=' ')[0:2].astype(int)

    # Read the sparse matrix with csr_matrix((data, (row_ind, col_ind)), shape=(M, N))
    # Subract -1 to rows and cols index because csr_matrix expects a 0 based index
    # Traspose counts matrix to have Cells as rows and Genes as cols as expected by AnnData objects
    spliced = sparse.csr_matrix((mtx[:,2], (mtx[:,0]-1, mtx[:,1]-1)), shape = shape).transpose()
    unspliced = sparse.csr_matrix((mtx[:,3], (mtx[:,0]-1, mtx[:,1]-1)), shape = shape).transpose()
    ambiguous = sparse.csr_matrix((mtx[:,4], (mtx[:,0]-1, mtx[:,1]-1)), shape = shape).transpose()
    
    # Load Genes and Cells identifiers
    obs = pd.read_csv(path+'/Velocyto/raw/barcodes.tsv',
                  header = None, index_col = 0)
    # Remove index column name to make it compliant with the anndata format
    obs.index.name = None
    
    var = pd.read_csv(path+'/Velocyto/raw/features.tsv', sep='\t',
                  names = ('gene_ids', 'feature_types'), index_col = 1)
    
    # Build AnnData object to be used with ScanPy and ScVelo
    adata = anndata.AnnData(X = X, obs = obs, var = var,
                        layers = {'spliced': spliced, 'unspliced': unspliced, 'ambiguous': ambiguous})
    adata.var_names_make_unique()
    
    # Subset Cells based on STAR filters
    selected_barcodes = pd.read_csv(path+'/Gene/filtered/barcodes.tsv', header = None)
    adata = adata[selected_barcodes[0]]
    
    return adata.copy()
```

alternative version with only counts and not velocity related matrices

```python
def buildAnndataFromStar(path):
    """Generate an anndata object from the STAR aligner output folder
	without velocity related informations"""
	# Load Read Counts
    X = sc.read_mtx(path+'/Gene/filtered/matrix.mtx')
	# Transpose counts matrix to have Cells as rows and Genes as cols as expected by AnnData objects
    X = X.X.transpose()
	# Load Genes and Cells identifiers
    obs = pd.read_csv(path+'/Gene/filtered/barcodes.tsv',
                  header = None, index_col = 0)
    # Remove index column name to make it compliant with the anndata format
    obs.index.name = None
    
    var = pd.read_csv(path+'/Gene/filtered/features.tsv', sep='\t',
                  names = ('gene_ids', 'feature_types'), index_col = 1)

	# Build AnnData object to be used with ScanPy
    adata = anndata.AnnData(X = X, obs = obs, var = var)
    adata.var_names_make_unique()
	return adata.copy()
```



**Load data** 

```python
NCAMpos = buildAnndataFromStar('data/prostate_scRNAseq/LNCaP_NCAM1.Solo.out/')
NCAMpos.obs['Origin'] = 'NCAM1+'

NCAMneg = buildAnndataFromStar('data/prostate_scRNAseq/LNCaP_NCAM1neg.Solo.out/')
NCAMneg.obs['Origin'] = 'NCAM1-'

NCAMbulk = buildAnndataFromStar('data/prostate_scRNAseq/LNCaP_bulk.Solo.out/')
NCAMbulk.obs['Origin'] = 'bulk'
```

**Concatenate datasets**

```python
adata = NCAMpos.concatenate(NCAMneg, NCAMbulk)
```

**Export merged dataset**

So you can restart from this step for re-runs

```python
adata.write('processed/LNCaP_merged_raw.h5ad', compression = 'gzip')
adata = anndata.read_h5ad('processed/LNCaP_merged_raw.h5ad')
```

## Quality check

```python	
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

adata.var['ribosomal'] = adata.var_names.str.contains('^RP[SL]')  # annotate the group of ribosomal RNA
sc.pp.calculate_qc_metrics(adata, qc_vars=['ribosomal'], percent_top=None, log1p=False, inplace=True)
```

Plot quality metrics by different merged experiments

```python
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True, groupby='Origin')
plt.hist(adata.obs['pct_counts_mt'])
plt.hist(adata.obs['pct_counts_ribosomal'])
#plt.hist(adata.obs['pct_counts_mt'])
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
sc.pl.scatter(adata, x='pct_counts_mt', y='n_genes_by_counts')
```

## Data preprocessing

**Subset cells**

```python
adata = adata[adata.obs['pct_counts_mt'] <= 20]
adata = adata[adata.obs['n_genes_by_counts'] >= 500]
adata = adata[adata.obs['total_counts'] >= 2000]
# Optionaly filter out cells with an average number of reads for each gene >= 6
# adata = adata[(adata.obs['total_counts']/adata.obs['n_genes_by_counts']) >= 6]

print('Number of cells passing filter: ' + str(adata.n_obs))
```

**Remove mitochondrial prior to normalization**

```python
# Remove reads mapped to mitochondrial DNA in order to scale after mitochondrial DNA removal
adata = adata[:, adata.var['mt'] == False]
```



## Analyze PCA Loadings

Pc loadings should be computed on the whole list of genes and not only highly variable ones `sc.pp.pca(adata, n_comps=50, use_highly_variable=True)` to analyze pathways

```python
PC_loading = pd.DataFrame(index = adata.var_names)
for PC in range(0, adata.varm['PCs'].shape[1]):
	PC_loading['PC' + str(PC+1)] = adata.varm['PCs'][:,PC]
```

Analyze with GSEA

```python
import gseapy as gp
from gseapy.plot import barplot, dotplot
# gp.get_library_name()

selected_PC = 'PC1'
pre_res = gp.prerank(rnk=PC_loading.loc[:,['symbol', selected_PC]].sort_values(selected_PC, ascending=False),
					 gene_sets='MSigDB_Hallmark_2020', #KEGG_2021_Human
                     processes=8,
                     permutation_num=100,
                     seed=6,
                     no_plot=True)

pre_res.res2d.query('fdr < 0.05').sort_values('nes', ascending=False)
```

### Analyze with ORA

First explore the distribution of loadings to select best elbow cutoff for number of genes to analyse

```python
selected_PC = 'PC1'
plt.stem(list(range(500)), PC_loading.sort_values(selected_PC, ascending = False)[selected_PC][:500])
plt.show()
plt.stem(list(range(500)), PC_loading.sort_values(selected_PC, ascending = False)[selected_PC][-500:])
plt.show()
```

Run PCA

```python
selected_PC = 'PC1'
enr = gp.enrichr(gene_list= PC_loading.sort_values(selected_PC, ascending=False)['symbol'][:50].to_list(),
                 gene_sets=['MSigDB_Hallmark_2020'], #'KEGG_2021_Human', 
                 organism='Human', 
                 description='test_name',
                 no_plot = True,
                 cutoff=0.05 # test dataset, use lower value from range(0,1)
                )
dotplot(enr.res2d)
```

## Plot PCs

```python
PC_projections = pd.DataFrame(index = adata.obs_names)
for PC in range(0, adata.obsm['X_pca'].shape[1]):
	PC_projections['PC' + str(PC+1) + '_position'] = adata.obsm['X_pca'][:,PC]

adata.obs = pd.concat([adata.obs, PC_projections], axis = 1)
```



## run gsea of marker genes for each cluster

```python
def extract_markers(adata):
	# initialize empty result
	result = pd.DataFrame({'Gene': [],
						   'Cluster': [],
						   'log2FC': [],
						   'pval': []})

	for cluster in adata.uns['rank_genes_clusters']['logfoldchanges'].dtype.names:
		partial_result = pd.DataFrame({'Gene': adata.uns['rank_genes_clusters']['names'][cluster],
				'Cluster': cluster,
				'log2FC': adata.uns['rank_genes_clusters']['logfoldchanges'][cluster],
				'pval': adata.uns['rank_genes_clusters']['pvals'][cluster],
				'pvals_adj': adata.uns['rank_genes_clusters']['pvals_adj'][cluster]})
		result = pd.concat([partial_result, result], axis=0)
	
	# result = result.query('pvals_adj < 0.05')
		
	return result 
```

```python
# run the analysis
sc.tl.rank_genes_groups(adata, groupby='leiden', key_added='rank_genes_clusters', method='wilcoxon', use_raw=False)
marker_genes = extract_markers(adata)
marker_genes['score'] = marker_genes['log2FC'] * (- np.log10(marker_genes['pval']))

# Export
with pd.ExcelWriter("result/AIRC-report/NCAM1pos-DE_genes/LNCaP_NCAMpos-Marker_gene_by_cluster.xlsx") as writer:
    marker_genes.to_excel(writer, sheet_name="marker_genes", index=False, freeze_panes = (True, False)) 
```



```python
import gseapy as gp
from gseapy.plot import barplot, dotplot

# initialize empty dataframe
result = pd.DataFrame({'Cluster': [], 'Term': [], 'ES': [], 'NES': [], 'NOM p-val': [], 'FDR q-val': [], 'Lead_genes': []})

# iterate for each cluster
for cluster in adata.uns['rank_genes_clusters']['logfoldchanges'].dtype.names:

	# run gsea preranked for each fluster sorting by log2FC
	prerank_res = gp.prerank(rnk=marker_genes.query('Cluster == "' + cluster + '"').loc[:,['Gene','log2FC']].sort_values('log2FC', ascending=False),
							 gene_sets='KEGG_2021_Human', #'MSigDB_Hallmark_2020'
                     		 processes=8,
                     		 permutation_num=100,
                     		 seed=6,
                     		 no_plot=True)
	# filter significant pathway and keep relevant column only
	partial_res = prerank_res.res2d.query('`FDR q-val` < 0.05').sort_values('NES', ascending=False).loc[:,['Term', 'ES', 'NES', 'NOM p-val', 'FDR q-val', 'Lead_genes']]
	
	# add a column specifying the cluster
	partial_res['Cluster'] = cluster

	result = pd.concat([partial_res, result], axis=0)
```

