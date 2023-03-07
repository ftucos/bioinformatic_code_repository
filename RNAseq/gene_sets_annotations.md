# Gene Sets Annotations

Useful annotated genesets for GSEA/ORA

### MSigDB

Main resource used for every annotation

https://www.gsea-msigdb.org/gsea/msigdb/

To retreive the annotation set within R you can use:

```R
library(msigdbr)
signature = msigdbr(species = "human")
```

or 

```R
library(msigdf)
msigdf::msigdf.human
# or 
msigdf::msigdf.mouse
```

there is some discrepancy within the two libraries, this is probabile due to the presence in the msigdf retreived annotation sets of the hortologs of mouse annotation sets

**NOTE:** you should not use Msigdb for kegg annotations since they are not updated since 2011 and pathways annotations doubled since then. You can use the `enrichKEGG()` funcation that queries latest annotations or manually load latest kegg annotations

## Manually retreive latest kegg annotation

In the following example we will be using the rest api to retreive the mapping withing each pathway and gene names (including aliases)

```r
# extract valid gene names
ensembl2symbol <- biomaRt::getBM(mart = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl"),
                                 attributes = c("external_gene_name", "external_synonym", "ensembl_gene_id", "gene_biotype", "description")) %>%
  ## Remove source info in the description field
  mutate(description = str_remove(description, "\\s?\\[.*\\]\\s?"))

alias2hcgn.raw <- ensembl2symbol %>%
  select(external_gene_name, external_synonym) %>%
  filter(!is.na(external_gene_name) & !is.na(external_synonym)) %>%
  distinct()

alias2hcgn <- set_names(x=alias2hcgn.raw$external_gene_name, nm=alias2hcgn.raw$external_synonym )
```



```R
library(tidyverse)

pathway2gene = read_tsv("http://rest.kegg.jp/link/mmu/pathway", col_names = c("pathway_id", "gene_id"))
pathway2name = read_tsv("http://rest.kegg.jp/list/pathway/mmu", col_names = c("pathway_id", "pathway"))
gene2name = read_tsv("http://rest.kegg.jp/list/mmu", col_names = c("gene_id", "gene")) 

gene2name <- gene2name %>%
  # gene with no semycolons have no name (eg. annotated cDNA sequences)
  filter(str_detect(gene, ";")) %>%
  mutate(gene = str_extract(gene, "^[^;]+") %>% strsplit(",\\s")) %>%
  # multiple gene names mapping the same id beacuse the mapping is based around orthologues
  unnest(cols = "gene") %>%
  # update alias names
  mutate(external_gene_name = case_when(
    gene_name %in% ensembl2symbol$external_gene_name ~ gene_name,
    gene_name %in% ensembl2symbol$external_synonym ~ alias2hcgn[gene_name],
    TRUE ~ as.character(NA)
  )) %>%
  select(-gene_name) %>%
  filter(!is.na(external_gene_name)) %>%
  distinct()

# no empty fields
table(is.na(gene2name))

kegg_t2g <- pathway2name %>%
  left_join(pathway2gene) %>%
  left_join(gene2name)

table(is.na(kegg_t2m))

# remove 11 pathways with missing genes
kegg_t2g <- kegg_t2g %>%
  filter(!is.na(gene))

table(is.na(kegg_t2m))

kegg_t2g <- kegg_t2g %>%
  select(gene, pathway) %>%
  mutate(pathway = str_remove(pathway, " - Mus musculus \\(mouse\\)"))

write.table(kegg_t2g, "processed/KEGG_mmu_term2gene.tsv", sep = "\t", row.names = F)

```



### Mitocarta for metabolic annotation

https://www.broadinstitute.org/mitocarta/mitocarta30-inventory-mammalian-mitochondrial-proteins-and-pathways

### YAP-TAZ signature

Hippo/YAP-TAZ pathway is difficoult to retreive from common GSEA (even absent in most annotation sets) because it's regulation is mainly post-translational. 

A custom signature of specific downstream targets has been defined:

> Wang, Yumeng et al. “Comprehensive Molecular Characterization of the Hippo Signaling Pathway in Cancer.” Cell reports vol. 25,5 (2018): 1304-1317.e5. doi:10.1016/j.celrep.2018.10.001

```R
YAPTAZ_targets_signature <- c("MYOF", "AMOTL2", "LATS2", "CTGF", "CYR61", "ANKRD1", "ASAP1", "AXL", "F3", "IGFBP3", "CRIM1", "FJX1", "FOXF2", "GADD45A", "CCDC80", "NT5E", "DOCK5", "PTPN14", "ARHGEF17", "NUAK2", "TGFB2", "RBMS3")
```

### Collection of pathway Gene Sets

http://baderlab.org/GeneSets

They keep updating genesets periodically from different resources, offering also an "Human_allpathways*.gmt" with all the resources grouped together. This resource is particularry convenient for KEGG annotation. Some years ago KEGG changed the redistribution policy for it's gene sets and now manual download from theri website si required. Barderlab does that (don't think it's legit).

## Manually parse GMT files

```R
parse_gmt_list <- function(element, gmt) {
  data.frame(gs_name = element,
             external_gene_name = gmt[[element]])
}

import_gmt <- function(path) {
    gmt <-  qusage::read.gmt(path)
    gmt_elements <- names(gmt)
  
   result <- map(gmt_elements, ~parse_gmt_list(., gmt)) %>%
     do.call(what=rbind) %>%
     left_join(ensembl2symbol %>% select(external_gene_name, ensembl_gene_id)) %>%
     select(-external_gene_name, ensembl_gene = ensembl_gene_id) %>%
     filter(!is.na(ensembl_gene)) %>%
     distinct() 
   result
}

pathway <- import_gmt("HumanCyc_2016.gmt")
```

