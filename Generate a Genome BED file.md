# Generate a Genome BED file

### Gene BED file from GTF

`makeTxDbFromGFF` supports also .gz comrpessed gtf files

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

#### with BedOps

1. the first awk skips comment lines
2. the second extracts gene levels annotations only
3. the two sed reformats the name field to avoid spaces and fix ambiguities betwen tab separator that you will encounter if using bedops tools

```bash	
awk '/^[^#]/ { print $0 }' gencode.v38.primary_assembly.annotation.gtf | awk 'BEGIN{FS="\t";OFS="\t"} $3 == "gene" { print $1,$4,$5,$9,$6,$7 }' | sed -e 's/; /;/g' | sed -e 's/ /=/g' > gencode.v38.primary_assembly.annotation.bed
```

use the generated file to annotate a bed file with the 'name' field (column 4)

```bash
bedmap --echo --echo-map-id --delim ';' --multidelim  '|' ranges.bed gencode.v38.primary_assembly.annotation.bed > ranges.annotated.bed
```

for control freec CNVs

```bash
bedmap --echo --echo-map-id --delim '|' --multidelim  '|' <(awk -F '\t' '{ print "chr"$1, $2, $3, $4 }' sample.tumor.mpileup.gz_CNVs) gencode.v38.primary_assembly.annotation.bed > sample.annotated.bed
```



