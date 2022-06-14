# Generate a Genome BED file

### Gene BED file from GTF

```R
library(AnnotationDbi)
library(GenomicFeatures)

txdb.GRCh38.v38 <- makeTxDbFromGFF("/Volumes/TucosHDD/Bioinformatics/resources/genomes/homosapiens/GENCODE/GRCh38/gencode.v38.primary_assembly.annotation.gtf",
                        format="gtf" , organism='Homo sapiens', dataSource="GRCh38 GENCODE v38 - primary_assembly")

saveDb(txdb.GRCh38.v38, file = "/Volumes/TucosHDD/Bioinformatics/resources/genomes/homosapiens/GENCODE/GRCh38/txdb.GRCh38.v38.sqlite")
txdb.GRCh38.v38 <- loadDb(file = "/Volumes/TucosHDD/Bioinformatics/resources/genomes/homosapiens/GENCODE/GRCh38/txdb.GRCh38.v38.sqlite")

genes <- genes(txdb.GRCh38.v38)
bed <- as.data.frame(genes) %>%
    # keep only chr1-22 and chrX genes
    filter(str_detect(seqnames, "chr[0-9X]{1,2}")) %>%
    transmute(chrom = seqnames, chromStart = start, chromEnd = end, name = gene_id, score = 0, strand = strand)

write_tsv(bed, "processed/genome/GRCh38.v38_genes.bed", col_names = F)
```

#### From UCSC

```R
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
hg38 <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes <- genes(hg38)
# ... follow previous instructions
```

