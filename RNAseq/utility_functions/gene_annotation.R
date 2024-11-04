##############################
#           OPTIONS          #
##############################

useCache <- TRUE
cacheDir <- "/Users/tucos/.bioinfoCache"

# pick the latest chache files available
biomartCache.human <- list.files(cacheDir, "ensembl2symbol-hsapiens_gene_ensembl.tsv", full.names = T) %>% sort(decreasing = T) %>% head(n=1)
biomartCache.mouse <- list.files(cacheDir, "ensembl2symbol-mmusculus_gene_ensembl.tsv", full.names = T) %>% sort(decreasing = T) %>% head(n=1)
t2gCache.human <- list.files(cacheDir, "t2g-human.tsv.gz", full.names = T) %>% sort(decreasing = T) %>% head(n=1)
t2gCache.mouse <- list.files(cacheDir, "t2g-mouse.tsv.gz", full.names = T) %>% sort(decreasing = T) %>% head(n=1)

date <- Sys.Date()

##############################
#             RUN            #
##############################

library(data.table)
library(tidyverse)
library(glue)
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("rename", "dplyr")

if (species %in% c("human", "homo sapiens", "hsapiens", "hsa")) {
  biomartCache <- biomartCache.human
  t2gCache <- t2gCache.human
  
  biomart_species = "hsapiens_gene_ensembl"
  kegg_species = "hsa"
  
} else if (species %in% c("mouse", "mus musculus", "mmusculus", "mmu")) {
  biomartCache <- biomartCache.mouse
  t2gCache <- t2gCache.mouse
  
  biomart_species = "mmusculus_gene_ensembl"
  kegg_species = "mmu"
} else {
  stop(glue("Error: {species} not allowed. Species must be either 'human' or 'mouse'"))
  
}

# Gene ID conversion table -----------------------
if(useCache) {
  ensembl2symbol = fread(biomartCache)
} else {
  mart <- biomaRt::useMart("ensembl", dataset = biomart_species)
  ensembl2symbol <- biomaRt::getBM(mart = mart,
                                   attributes = c("external_gene_name", "ensembl_gene_id",
                                                  "gene_biotype",
                                                  "description")) %>%
    ## Remove source info in the description field
    mutate(description = str_remove(description, "\\s?\\[.*\\]\\s?"))
  write_tsv(ensembl2symbol, paste0(cacheDir, "/", date, "-ensembl2symbol", "-", biomart_species, ".tsv"))
}

# GSEA annotation data -----------------------------------------------------------------
# Load latest KEGG pathway annotation data -----------------------------------
if(useCache) {
  t2g = fread(t2gCache)
} else {
  pathway2gene = read_tsv(paste0("http://rest.kegg.jp/link/",kegg_species,"/pathway"), col_names = c("pathway_id", "gene_id"))
  pathway2name = read_tsv(paste0("http://rest.kegg.jp/list/pathway/", kegg_species), col_names = c("pathway_id", "pathway"))
  gene2name = read_tsv(paste0("http://rest.kegg.jp/list/", kegg_species), col_names = c("gene_id", "biotype", "X3", "gene")) 
  
  gene2name <- gene2name %>%
    # gene with no semicolons have no name (eg. annotated cDNA sequences)
    filter(str_detect(gene, ";")) %>%
    mutate(gene = str_extract(gene, "^[^;]+") %>% strsplit(",\\s")) %>%
    # multiple gene names mapping the same id beacuse the mapping is based around orthologues
    unnest(cols = "gene")
  
  # no empty fields
  #table(is.na(gene2name))
  
  kegg_t2g <- pathway2name %>%
    left_join(pathway2gene %>% mutate(pathway_id = str_remove(pathway_id, "path:"))) %>%
    left_join(gene2name)
  
  table(is.na(kegg_t2g))
  
  # remove 11 pathways with missing genes
  kegg_t2g <- kegg_t2g %>%
    filter(!is.na(gene))
  
  #table(is.na(kegg_t2g))
  
  table(kegg_t2g$gene %in% ensembl2symbol$external_gene_name)
  
  kegg_t2g <- kegg_t2g %>%
    mutate(pathway = str_remove(pathway, " - Homo sapiens \\(human\\)") %>% str_remove(" - Mus musculus \\(house mouse\\)"))%>%
    # adapt to the msigdb column names
    select(external_gene_name = gene, gs_name = pathway) %>%
    inner_join(ensembl2symbol %>% select(external_gene_name, ensembl_gene = ensembl_gene_id), by="external_gene_name") %>%
    distinct() %>%
    mutate(gs_subcat = "KEGG_latest") %>%
    dplyr::rename("gene_symbol" = "external_gene_name")
  
  ## Load annotation sets for GSEA
  msigdb <- msigdbr(species = species)
  # clean version of the annotation database
  t2g <- msigdb %>%
    mutate(gs_subcat = ifelse(gs_cat == "H", "HALLMARK", str_remove(gs_subcat, "^CP:"))) %>%
    # select signature of interest
    filter(gs_subcat %in% c("HALLMARK", "KEGG", "REACTOME", "PID", "BP", "WIKIPATHWAYS", "GO:BP")) %>%
    mutate(gs_name = case_when(
      gs_subcat == "HALLMARK" ~ str_remove(gs_name, "^HALLMARK_") %>% str_replace_all("_", " "),
      gs_subcat == "KEGG" ~ gs_description,
      gs_subcat == "REACTOME" ~ gs_description,
      gs_subcat == "PID" ~ gs_description,
      gs_subcat == "WIKIPATHWAYS" ~ gs_description,
      gs_subcat == "GO:BP" ~ str_remove(gs_name, "^GOBP_") %>% str_replace_all("_", " "),
      TRUE ~ gs_name)) %>%
    # add latest kegg annotations
    bind_rows(kegg_t2g)
  
  write_tsv(t2g, paste0(cacheDir, "/", date, "-t2g", "-", species, ".tsv"))
  # compress the file
  system(paste0("gzip ", cacheDir, "/", date, "-t2g", "-", species, ".tsv"))
  
  # remove intermediate files 
  rm(gene2name, kegg_t2g, pathway2gene, pathway2name, biomart_species, kegg_species, date)
}

# Pathway annotation table --------------------------------------------
pathway_annotation_table <- full_join(
  t2g %>% filter(gs_subcat == "HALLMARK") %>%
    select(gs_name, ensembl_gene) %>%
    group_by(ensembl_gene) %>%
    summarize(pathways = paste0(gs_name, collapse = "|")) %>%
    dplyr::rename("ensembl_gene_id" = "ensembl_gene",
                  "HALLMARK" = "pathways"),
  t2g %>% filter(gs_subcat == "KEGG_latest") %>%
    group_by(ensembl_gene) %>%
    summarize(pathways = paste0(gs_name, collapse = "|")) %>%
    dplyr::rename("ensembl_gene_id" = "ensembl_gene",
                  "KEGG_latest" = "pathways")) %>%
  full_join(
    t2g %>% filter(gs_subcat == "REACTOME") %>%
      group_by(ensembl_gene) %>%
      summarize(pathways = paste0(gs_name, collapse = "|")) %>%
      dplyr::rename("ensembl_gene_id" = "ensembl_gene",
                    "REACTOME" = "pathways")) %>%
  full_join(
    t2g %>% filter(gs_subcat == "PID") %>%
      group_by(ensembl_gene) %>%
      summarize(pathways = paste0(gs_name, collapse = "|")) %>%
      dplyr::rename("ensembl_gene_id" = "ensembl_gene",
                    "PID" = "pathways"))


