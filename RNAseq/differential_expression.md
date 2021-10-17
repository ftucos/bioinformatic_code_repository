# Differential Expression

## Import counts with TxImport

To import counts from TxImport from transcript quantifier

### Salmon

```R
library(tximport)
library(GenomicFeatures)
txdb <- makeTxDbFromGFF("/Volumes/TucosHDD/Bioinformatics/resources/genomes/homosapiens/GENCODE/GRCh38/gencode.v38.primary_assembly.annotation.gtf",
                                            format="gtf" , organism='Homo sapiens', dataSource="GRCh38 GENCODE v38 - primary_assembly")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME") 
```

```R
# Import with TxImport ---------------------------------------
filedir <- "/Volumes/TucosHDD/Bioinformatics/data/breast_CDK12_OE/BulkRNAseq-CDK12_MCF10A/salmon"
files <- list.files(filedir, pattern="*quant.sf", recursive = T, full.names = T)
names(files) <- str_extract(files, "(?<=salmon\\/)[^/]+")

## Import gene abundances and converts them to counts (integers)
txi <- tximport(files, type = "salmon", txIn = T, txOut = F,  tx2gene = tx2gene)
names(txi)

table(txi$length != 0)

## Remove version from ensembl gene id
rownames(txi$abundance) <- str_remove(rownames(txi$abundance), "\\..*")
rownames(txi$counts) <- str_remove(rownames(txi$counts), "\\..*")
rownames(txi$length) <- str_remove(rownames(txi$counts), "\\..*")
```

### RSEM

```R
filedir <- "/Volumes/TucosHDD/Bioinformatics/workspace/siNUMB_RNAseq_processing/results/star_rsem"
files <- list.files(filedir, pattern="*genes.results", recursive = T, full.names = T)
names(files) <- str_extract(files, "[^/]+(?=\\.genes\\.results$)")

## Import gene abundances and converts them to counts (integers)
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

## fix "Error: all(lengths > 0) is not TRUE" error
txi.rsem$length[txi.rsem$length == 0] = 0.01

## Remove version from ensembl gene id
rownames(txi.rsem$abundance) <- str_remove(rownames(txi.rsem$abundance), "\\..*")
rownames(txi.rsem$counts) <- str_remove(rownames(txi.rsem$counts), "\\..*")
rownames(txi.rsem$length) <- str_remove(rownames(txi.rsem$counts), "\\..*")

## Remove "-" character from col names

colnames(txi.rsem$abundance) <- str_remove(colnames(txi.rsem$abundance), "-")
colnames(txi.rsem$counts) <- str_remove(colnames(txi.rsem$counts), "-")
colnames(txi.rsem$length) <- str_remove(colnames(txi.rsem$counts), "-")
```



## edgeR Differential Expression

```R
DEG <-  DGEList(txi$counts, group=metadata$condition)
keep <- filterByExpr(DEG) 
DEG <- DEG[keep, ,keep.lib.sizes=T]
DEG <- calcNormFactors(DEG)
plotMDS(DEG)
design <- model.matrix(~metadata$condition + metadata$replicate )
DEG <- estimateDisp(DEG, design) # robust = T

fit <- glmQLFit(DEG, design) #  robust = T
# coeff 2 = metadata$conditionCDK12
qlf <- glmQLFTest(fit, coef=2)
res <- topTags(qlf, n = "Inf")
```

Use `tr <- glmTreat(fit, coef=2, lfc=0.5)` in place of `glmQLFTest()` in order to compute p-values for the difference being higher than the logFC specified rather than just different from 0.



```R
result <- res$table %>%
  rownames_to_column("ensembl_gene_id") %>%
  left_join(ensembl2symbol) %>%
  # rename columns to match them with the ones obtained from DEseq2 to reuse the same function
  select(ensembl_gene_id, symbol = external_gene_name, log2FoldChange = logFC, pvalue = PValue, padj = FDR, mean_gene_abundance = logCPM, entrezgene_id, gene_biotype, description) %>%
  mutate(DEscore = log2FoldChange * -log10(pvalue)) %>%
  mutate(padj = p.adjust(pvalue,method = "hommel"))
```

## DEseq2

```R
dds <- DESeqDataSetFromTximport(txi = txi,
                                colData = metadata,
                                design = ~ condition + replicate) # batch correction
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="condition_CDK12_vs_EV")
# shrinked results for low counts 
res.shrink <- lfcShrink(dds, coef="condition_CDK12_vs_EV", type="apeglm", svalue = FALSE)
```

You can pass the `lfcThreshold = 1` argument to `result()` in order to compute p-values for the difference being higher than the logFC specified rather than just different from 0.

### Extract counts

```R
dds <- estimateSizeFactors(dds)
counts <- counts(dds, normalized=T)
```

### Remove batch effect

To be performed on vst or rlog transformed data or logCPM in case of edgeR

```R
dds <- DESeqDataSetFromMatrix(countData=counts, colData=factors, design = ~ Batch + Covariate)
dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE)
mat <- assay(vsd)
mat <- limma::removeBatchEffect(mat, vsd$Batch)
assay(vsd) <- mat
counts_batch_corrected <- assay(vsd)
```

### Custom PCA function for DEseq2

```r 
library(ggrepel)

# calculate the variance for each gene
rv <- rowVars(assay(vsd))

# select the ntop genes by variance
# you can select them just by ranking by variance since a variance stabilizing tra nsformation
# has been aleready applied
select <- order(rv, decreasing=TRUE)[seq_len(min(1000, length(rv)))]

PCA.raw <- prcomp(t(assay(vsd[select,])), center=T, scale = T)
PCA <- PCA.raw$x %>%
  as.data.frame() %>%
  select(PC1, PC2, PC3, PC4, PC5) %>%
  rownames_to_column("sample") %>%
  left_join(colData)

# the contribution to the total variance for each component
PCA.var <- PCA.raw$sdev^2 / sum(PCA.raw$sdev^2)

ggplot(PCA, aes(x=PC1, y=PC2, label=sample, color=pct.usable.reads))+
  geom_point(size=2)+
  geom_text_repel(max.overlaps = 30, min.segment.length = 0.1)+
  scale_color_continuous(high="green", low="red")+
  theme_bw()+
  xlab(paste0("PC1 (", round(PCA.var[1]*100, 1), "%)"))+
  ylab(paste0("PC2 (", round(PCA.var[3]*100, 1), "%)")) +
  ggtitle("Identification of samle outliers by RNA quality")
```

