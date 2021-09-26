# How to handle TCGA data portal

1. Select cases and files lists from GDC Data Portal https://portal.gdc.cancer.gov/
2. Add selected files to the chart
3. from the cart you can dowload the manifest, clinical data and samplesheet (matching file-id to case-id)

## Dowload TCGA data from manifest

```shell
#!/bin/bash
set ue -o pipefail

WD="/Volumes/TucosHDD/Bioinformatics/workspace/TCGA-bladder-NUMB_analysis"
cd $WD

mkdir -p $WD/data/RNAseq_counts
cd $WD/data/RNAseq_counts
gdc-client download -m  $WD/data/manifests/gdc_manifest.RNAseq_counts.txt

mkdir -p $WD/data/mirnaRNAseq_quant
cd $WD/data/RNAseq_counts
gdc-client download -m  $WD/data/manifests/gdc_manifest.miRNAseq_quantification.txt
```

# Merge reads into one table 

```r
setwd("/Volumes/TucosHDD/Bioinformatics/workspace/TCGA-bladder-NUMB_analysis")
library('data.table')
library('tidyverse')

# RNAseq unificataion -----------------------------------------------------------------------
# Samplesheet with filename to patient/sample-id matching
RNAseq_samplesheet <- read_delim("data/metadata/RNAseq/gdc_sample_sheet.2021-06-27.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
# Generate a named vector for easily converting filename to sample-id 
filename_to_sampleID <- set_names(RNAseq_samplesheet$`Sample ID`, RNAseq_samplesheet$`File Name`)

# Function for parsing read counts
parse_RNAseq_reads <- function(path) {
  file_name <- str_extract(path, "(?<=/)[^/]*$")
  counts <- fread(path, col.names = c("ensembl_gene_id", filename_to_sampleID[file_name])) %>%
    column_to_rownames("ensembl_gene_id")
}

# List files to parse
RNAseq_files <- list.files("data/RNAseq_counts", pattern = ".htseq.counts.gz$", recursive = T, full.names = T)

# parse all read counts and join them in a single table
RNAseq <- map(RNAseq_files, parse_RNAseq_reads) %>%
  do.call(what = cbind)

# remove entire sample with duplicated sample-id
RNA_duplicated_samples <- colnames(RNAseq)[colnames(RNAseq) %>% duplicated()]
RNAseq <- RNAseq[,!(colnames(RNAseq) %in% RNA_duplicated_samples)]

# remove version from ensembl_gene_id
RNAseq.1 <- RNAseq %>% 
  rownames_to_column("ensembl_gene_id") %>%
  mutate(ensembl_gene_id = str_remove(ensembl_gene_id, "\\..*$"))

# check no duplicated gene id are present after version removal
RNAseq.1$ensembl_gene_id %>% duplicated() %>% table()

write_tsv(RNAseq.1, file = "processed/TCGA-BLCA_RNAseq_counts.tsv",quote_escape = F)


# miRNAseq unificataion -----------------------------------------------------------------------
miRNAseq_samplesheet <- read_delim("data/metadata/miRNAseq/gdc_sample_sheet.2021-06-27.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
filename_to_sampleID <- set_names(miRNAseq_samplesheet$`Sample ID`, miRNAseq_samplesheet$`File Name`)


# Function for parsing read counts
parse_miRNA_reads <- function(path) {
  file_name <- str_extract(path, "(?<=/)[^/]*$")
  counts <- fread(path, col.names = c("miRNA", filename_to_sampleID[file_name], "scaled", "cross-mapped")) %>%
    select(-scaled, -`cross-mapped`) %>%
    column_to_rownames("miRNA")
}

# List files to parse
miRNAseq_files <- list.files("data/RNAseq_counts", pattern = "mirnas.quantification.txt$", recursive = T, full.names = T)

# parse all read counts and join them in a single table
miRNAseq <- map(miRNAseq_files, parse_miRNA_reads) %>%
  do.call(what = cbind)

# remove entire sample with duplicated sample-id
miRNA_duplicated_samples <- colnames(miRNAseq)[colnames(miRNAseq) %>% duplicated()]
miRNAseq <- miRNAseq[,!(colnames(miRNAseq) %in% miRNA_duplicated_samples)] %>%
  rownames_to_column("miRNA")

write_tsv(miRNAseq, file = "processed/TCGA-BLCA_miRNAseq_counts.tsv",quote_escape = F)
```

**Note:** In the following code I'll remove any duplicated sample id since I want to correlate transcripts to miRNAs and I prefer to be on the safe side in selecting the proper cases/aliquots. For simple transcriptomics analysis, you could rather sum the reads or keep the one with the highest non 0 counts number of genes.
