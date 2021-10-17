# Genomic Annotation

Note: most of those libraries will overwrite the dplyr `select()` function, so you have to load dplyr as last library.

## Gene level annotations

### biomaRt

Build a remapping data frame

```R
ensembl2symbol <- biomaRt::getBM(values = c("ENSG00000012048", "ENSG00000111328"),
                                 filters = "ensembl_gene_id",
                                 mart = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl"),
                                 attributes = c("external_gene_name", "entrezgene_id", "ensembl_gene_id", "gene_biotype", "description"),
                                ) %>%
  ## Remove source info in the description field
  mutate(description = str_remove(description, "\\s?\\[.*\\]\\s?"))

```

If you want to gather the annotations available for every gene you can simply remove the `values` and `filters` arguments.

To list all the available annotations use `biomaRt::listAttributes(mart)`

```R
> head(biomaRt::listAttributes(mart))
                           name                  description         page
1               ensembl_gene_id               Gene stable ID feature_page
2       ensembl_gene_id_version       Gene stable ID version feature_page
3         ensembl_transcript_id         Transcript stable ID feature_page
```

### org.* packages

For gene level annotation

```R
library(org.Hs.eg.db)
ensembl2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                        keys =  c("ENSG00000012048", "ENSG00000111328"),
                                        columns = c("SYMBOL", "GENETYPE"),
                                        keytype = "ENSEMBL")
```

Or you can use it to return the matching vector (eg, within a `mutate()` function)

```R
genes <- data.frame(ensembl_gene_id = c("ENSG00000012048", "ENSG00000111328"))
genes <- genes %>%
  mutate(symbol = AnnotationDbi::mapIds(org.Hs.eg.db,
                                        keys = ensembl_gene_id,
                                        column = "SYMBOL", keytype = "ENSEMBL"))
```

To list all the available annotations use `columns(org.Hs.eg.db)`

```R
 [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"    
 [7] "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"     "GENETYPE"     "GO"          
[13] "GOALL"        "IPI"          "MAP"          "OMIM"         "ONTOLOGY"     "ONTOLOGYALL" 
[19] "PATH"         "PFAM"         "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"      
[25] "UCSCKG"       "UNIPROT" 
```

## Mouse/Human orthologous interconversion

```R
human <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
mouse2human <- biomaRt::getLDS(attributes = c("mgi_symbol"), mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
human2mouse <- biomaRt::getLDS(attributes = c("hgnc_symbol"), mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
```

Adapted from: https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/

## Transcript level annotations

### TxDb.* packages

For transcript level annotation and genomic features

It supports the same syntax of org.* packages since it uses AnnotationDbi under the hood

Example use case for plotting the NUMB transcripts

```R
library(org.Hs.eg.db) # To map gene symbol to entrezid
library(TxDb.Hsapiens.UCSC.hg19.knownGene) # to extract transcript level informations
library(Gviz) # to plot transcripts

entrez_id <- mapIds(org.Hs.eg.db, "CDKN2A", "ENTREZID",
              "SYMBOL")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txid <- select(txdb, entrez_id, "TXNAME", "GENEID")[["TXNAME"]]
cds <- cdsBy(txdb, by="tx", use.names=TRUE)
CDKN2Acds <- cds[names(cds) %in% txid]
tx <- rep(names(CDKN2Acds), lengths(CDKN2Acds))
id <- unlist(CDKN2Acds)$cds_id
grt <- GeneRegionTrack(CDKN2Acds, name="CDKN2A", id=tx,
                       gene="CDKN2A", feature=tx, transcript=tx, exon=id)
plotTracks(list(GenomeAxisTrack(), grt))

```

![image-20210706093029284](/Users/tucos/Library/Application Support/typora-user-images/image-20210706093029284.png)

To list all the available annotations use `columns(TxDb.Hsapiens.UCSC.hg38.knownGene)`

```R
 [1] "CDSCHROM"   "CDSEND"     "CDSID"      "CDSNAME"    "CDSPHASE"   "CDSSTART"   "CDSSTRAND" 
 [8] "EXONCHROM"  "EXONEND"    "EXONID"     "EXONNAME"   "EXONRANK"   "EXONSTART"  "EXONSTRAND"
[15] "GENEID"     "TXCHROM"    "TXEND"      "TXID"       "TXNAME"     "TXSTART"    "TXSTRAND"  
[22] "TXTYPE" 
```

## Genomic Sequeces

### BSgenome.* Packages

```R
library(BSgenome.Hsapiens.UCSC.hg38)
# following the code of the previous example
extractTranscriptSeqs(Hsapiens, CDKN2Acds)
```



