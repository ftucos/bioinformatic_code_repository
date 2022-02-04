library(msigdbr)
library(clusterProfiler)
library(patchwork)
library(ggrepel)
library(egg)
library(gridExtra)
library(tidyverse)

###############################
#          OPTIONS            #
###############################
species = "mouse"
if (species == "human") {
  biomart_species = "hsapiens_gene_ensembl"
  kegg_species = "hsa"
} else if (species == "mouse") {
  biomart_species = "mmusculus_gene_ensembl"
  kegg_species = "mmu"
    
}

# Gene ID conversion table -----------------------
mart <- biomaRt::useMart("ensembl", dataset = biomart_species)
ensembl2symbol <- biomaRt::getBM(mart = mart,
                                 attributes = c("external_gene_name", "ensembl_gene_id",
                                                "gene_biotype",
                                                "description")) %>%
  ## Remove source info in the description field
  mutate(description = str_remove(description, "\\s?\\[.*\\]\\s?"))

# Annotate Result Table ---------------------------
ResToTable <- function(res, package) {
   if (package == "DESeq2") {
     # res obtained from DESeq2::result(dds) or DESeq2::lfcShrink(dds) functions 
     result.table <- res %>% 
       as.data.frame() %>%
       rownames_to_column("ensembl_gene_id") %>%
       left_join(ensembl2symbol) %>%
       select(external_gene_name, log2FoldChange, pvalue, padj, mean_gene_abundance = baseMean, gene_biotype, description, ensembl_gene_id)
   
    } else if (package == "edgeR") {
     # res obtained from edgeR::topTags(qlf) function
     result.table <- res$table %>%
       rownames_to_column("ensembl_gene_id") %>%
       left_join(ensembl2symbol) %>%
       select(external_gene_name, log2FoldChange = logFC, pvalue = PValue,  padj = FDR, mean_gene_abundance = logCPM, gene_biotype, description, ensembl_gene_id)
    } else {
     warning("package argument must be either 'DESeq2' or 'edgeR'!")
     return(NULL)
   }
}

# Volcano Plot -------------------------------------
VolcanoPlot <- function(result, title = element_blank(), thrLog2FC = 1, thrPadj = 0.05, pctLabel = 0.005, protein_coding_label_only = FALSE) {
  
  # requires you to have applied ResToTable()
  x_limits <- c(-max(abs(result$log2FoldChange), na.rm=T), max(abs(result$log2FoldChange), na.rm=T))
  
  label_cutof <- result %>%
    # plot labels only for genes within the threshold defined
    filter(padj <= thrPadj, abs(log2FoldChange) >= thrLog2FC) %>%
    # rank priority of genes to plot based on combination of log2FC and padj
    # annotate protein coding only 
    filter((!protein_coding_label_only) | (ensembl_gene_id %in% (ensembl2symbol %>% filter(gene_biotype == "protein_coding") %>% pull("ensembl_gene_id")))) %>%
    mutate(DEscore = abs(log2FoldChange * -log10(pvalue))) %>%
    pull(DEscore) %>%
    # extract only the fraction of labels to plot
    quantile(1-pctLabel)
  
  df <- result %>%
    # Color code base on up/downregulation state defined by custom thresholds
    mutate(color = case_when(
                             log2FoldChange >= thrLog2FC & padj <= thrPadj ~ "up",
                             log2FoldChange <= -thrLog2FC & padj <= thrPadj ~ "down"),
           # labels to be shown need to be matching upregulated genes and be within the filtering threshold
          lab = ifelse(abs(log2FoldChange) >= thrLog2FC &
                         padj <= thrPadj &
                         # annotate protein coding only 
                         (!protein_coding_label_only) | (ensembl_gene_id %in% (ensembl2symbol %>% filter(gene_biotype == "protein_coding") %>% pull("ensembl_gene_id"))) &
                         abs(log2FoldChange * -log10(pvalue)) >= label_cutof, external_gene_name, NA))
  
  ggplot(df, aes(x=log2FoldChange, y=-log10(padj), color = color)) +
    geom_point(data=filter(df, is.na(color)), aes(x=log2FoldChange, y=-log10(padj)), color="gray", size=1) +
    geom_point(size = 1) +
    geom_vline(xintercept = c(-thrLog2FC, thrLog2FC), linetype = "dashed")+
    geom_hline(yintercept = -log10(thrPadj), linetype = "dashed")+
    geom_text_repel(aes(label = lab), size=3, color="black", segment.alpha = 0.7, max.overlaps = 40)+
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(size = 1),
          legend.position = "none") +
    xlim(x_limits) +
    scale_color_manual(values = c("dodgerblue", "brown1"), na.translate = F) +
    xlab(expression(log["2"](Fold~Change))) + ylab(expression(-log["10"](Adjusted~p~Value))) + 
    ggtitle(title)
}


# Load latest KEGG pathway annotation data -----------------------------------

pathway2gene = read_tsv(paste0("http://rest.kegg.jp/link/",kegg_species,"/pathway"), col_names = c("pathway_id", "gene_id"))
pathway2name = read_tsv(paste0("http://rest.kegg.jp/list/pathway/", kegg_species), col_names = c("pathway_id", "pathway"))
gene2name = read_tsv(paste0("http://rest.kegg.jp/list/", kegg_species), col_names = c("gene_id", "gene")) 

gene2name <- gene2name %>%
  # gene with no semicolons have no name (eg. annotated cDNA sequences)
  filter(str_detect(gene, ";")) %>%
  mutate(gene = str_extract(gene, "^[^;]+") %>% strsplit(",\\s")) %>%
  # multiple gene names mapping the same id beacuse the mapping is based around orthologues
  unnest(cols = "gene")

# no empty fields
#table(is.na(gene2name))

kegg_t2g <- pathway2name %>%
  left_join(pathway2gene) %>%
  left_join(gene2name)

table(is.na(kegg_t2g))

# remove 11 pathways with missing genes
kegg_t2g <- kegg_t2g %>%
  filter(!is.na(gene))

#table(is.na(kegg_t2g))

table(kegg_t2g$gene %in% ensembl2symbol$external_gene_name)

kegg_t2g <- kegg_t2g %>%
  mutate(pathway = str_remove(pathway, " - Mus musculus \\(mouse\\)" %>% str_remove(" - Homo sapiens (human)"))) %>%
  # adapt to the msigdb column names
  select(external_gene_name = gene, gs_name = pathway) %>%
  inner_join(ensembl2symbol %>% select(external_gene_name, ensembl_gene = ensembl_gene_id), by="external_gene_name") %>%
  distinct() %>%
  mutate(gs_subcat = "KEGG_latest") %>%
  rename("gene_symbol" = "external_gene_name")


# GSEA ----------------------------------------
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

## Prepare imput data for GSEA
PrepareForGSEA <- function(result, .na.rm = FALSE) {
  
  # chose whether remove genes with log2FC undefined (lowly expressed) or set theire log2FC to 0 since they are not activated by the treatment
  if (.na.rm == TRUE) {
    result <- result %>%
      filter(!is.na(log2FoldChange))
  } else {
    result <- result %>%
      mutate(log2FoldChange = ifelse(is.na(log2FoldChange), 0, log2FoldChange))
  }
  
  # Generate a sorted named vector with ensembl as name and log2FC as value
  .input <- setNames(result$log2FoldChange, result$ensembl_gene_id) %>%
    sort(decreasing = TRUE)
  
  return(.input)
}

## Plot GSEA table
plotGSEA <- function(.GSEA, title = "", cutoff = 0.05, subgroup = "all", fixed_dimensions=F) {
  
  if (subgroup == "pos") {
    .GSEA <- .GSEA %>% filter(NES > 0)
  } else if (subgroup == "neg") {
      .GSEA <- .GSEA %>% filter(NES < 0)
  }
  plot <- .GSEA %>%
    filter(qvalues <= cutoff) %>%
    mutate(Description = factor(Description, levels = 
                                  .GSEA %>% arrange(NES) %>% pull(Description))) %>%
    ggplot(aes(x=Description, y=NES, fill=qvalues)) +
    geom_col() + 
    theme_bw() +
    scale_fill_continuous(type = "gradient",
                          low = "#E74C3C", high = "#F1C40F",
                          space = "Lab", na.value = "grey50", guide = "colourbar")  +
    coord_flip() +
    ggtitle(paste0("GSEA: ", title)) +
    xlab("") +
    ylab("Normalized Enrichment Score") +
    theme(panel.grid = element_blank(),
          plot.title.position = "plot")
  
  # Fixed dimensions based on number of pathways
  if(fixed_dimensions) {
   plot <-  set_panel_size(plot,
      # define plot hieght based on the numer of pathways to be plotted
      height=unit((.GSEA %>% filter(qvalues <= cutoff ) %>% nrow())/1.5, "cm"),
      width=unit(8, "cm")
    )
    # plot the plot
    plot(plot)
    # return the plot (for export preservig dimensions)
    plot
  } else {
    # just plot the plot
    plot(plot)
    # and return it for export
    plot
  }
}

plotGSEAwide <- function(.GSEA, title = "", cutoff = 0.05, subgroup = "all") {
  
  if (subgroup == "pos") {
    .GSEA <- .GSEA %>% filter(NES > 0)
  } else if (subgroup == "neg") {
    .GSEA <- .GSEA %>% filter(NES < 0)
  }
  
  .GSEA %>%
    filter(qvalues <= cutoff) %>%
    mutate(Description = factor(Description, levels = 
                                  .GSEA %>% arrange(desc(NES)) %>% pull(Description))) %>%
    ggplot(aes(x=Description, y=NES, fill=qvalues)) +
    geom_col() + 
    theme_bw() +
    scale_fill_continuous(type = "gradient",
                          low = "#E74C3C", high = "#F1C40F",
                          space = "Lab", na.value = "grey50", guide = "colourbar")  +
    #   coord_flip() +
    ggtitle(paste0("GSEA: ", title)) +
    xlab("") +
    ylab("Normalized Enrichment Score") +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          plot.margin=unit(c(0.2,0.2,0.2,2),"cm"),
          plot.title.position = "plot")
}



## Run GSEA
runGSEA <- function(result, annotation, title = "", cutoff = 0.05, plot = FALSE, toDataFrame = TRUE, custom_annotation=F, ...) {
  if(custom_annotation) {
    # custom t2g object to be passed into annotatation argument
    selected_t2g <- annotation
  } else {
    selected_t2g <- # extract t2g of interest based on string passed into annotation
      t2g %>% filter(gs_subcat == annotation) %>% select(gs_name, ensembl_gene)
  }
    
  .em <- GSEA(PrepareForGSEA(result),
              TERM2GENE = selected_t2g,
              pAdjustMethod = "fdr",
              by = "fgsea",
              pvalueCutoff = 1,
              eps = 0,
              seed = 111 
              )
  .em.summary <- as.data.frame(.em)
  
  # If not specified, use the object name as plot title
  title <- ifelse(title == "", deparse(substitute(result)), title)
  # PLot
  if (plot == TRUE) {
    plot(plotGSEA(.em.summary, title = title, cutoff = cutoff, ...))
  }
  # return the summary table
  if(toDataFrame == TRUE) {
    return(.em.summary)
  } else {
    return(.em)
  }
}

# Pathway annotation table --------------------------------------------
pathway_annotation_table <- full_join(
    t2g %>% filter(gs_subcat == "HALLMARK") %>%
      select(gs_name, ensembl_gene) %>%
      group_by(ensembl_gene) %>%
      summarize(pathways = paste0(gs_name, collapse = "|")) %>%
      rename("ensembl_gene_id" = "ensembl_gene",
             "HALLMARK" = "pathways"),
    t2g %>% filter(gs_subcat == "KEGG_latest") %>%
      group_by(ensembl_gene) %>%
      summarize(pathways = paste0(gs_name, collapse = "|")) %>%
      rename("ensembl_gene_id" = "ensembl_gene",
             "KEGG_latest" = "pathways")) %>%
  full_join(
    t2g %>% filter(gs_subcat == "REACTOME") %>%
      group_by(ensembl_gene) %>%
      summarize(pathways = paste0(gs_name, collapse = "|")) %>%
      rename("ensembl_gene_id" = "ensembl_gene",
             "REACTOME" = "pathways")) %>%
  full_join(
    t2g %>% filter(gs_subcat == "PID") %>%
      group_by(ensembl_gene) %>%
      summarize(pathways = paste0(gs_name, collapse = "|")) %>%
      rename("ensembl_gene_id" = "ensembl_gene",
             "PID" = "pathways"))

# Define custom GSEA plot function ------------
custom_gseaplot <- function(res, pathway) {
  statistics <- res %>% as.data.frame() %>% filter(Description == pathway)
  gseaplot(res, geneSetID = pathway, by = "runningScore", title = pathway, color.line = "black") + 
    theme_minimal() + 
    ylab("Enrichment Score") +
    xlab("")+
    labs(caption = paste0("NES: ", round(statistics$NES,2),", p-value: ", round(statistics$pvalue,4), ", q-value: ",  round(statistics$qvalues,4)))
}

# Plot heatmap -----------------------
plotHEATMAP <- function(zscore, selected_genes, title, ylab) {
  heatmap.input <- zscore %>%
    filter(external_gene_name %in% selected_genes,
           !is.na(zscore))
  
  heatmap.input <- heatmap.input %>%
    mutate(external_gene_name = factor(external_gene_name,
                                       levels=result %>% arrange(log2FoldChange*-log10(pvalue)) %>% filter(external_gene_name %in% selected_genes) %>% pull("external_gene_name") %>% unique()))
  
  heatmap <- heatmap.input %>%
    mutate(sample = str_remove(sample, paste0(selected_sample, "_")) %>% str_replace_all("_", " ")) %>%
    ggplot(aes(x=sample, y=external_gene_name, fill=zscore)) +
    geom_tile()+
    scale_x_discrete(expand=c(0, 0)) + 
    scale_y_discrete(expand=c(0, 0))+
    scale_fill_distiller(type = "div",palette = 5, name="z-score(logCPM)",
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
    theme(legend.position = "none")+
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title.position = "plot",
          panel.border = element_rect(size=0.75),
          legend.position = "bottom")+
    ylab(ylab)+
    xlab("") +
    ggtitle(title)
  
  heatmap <-  set_panel_size(heatmap,
                             # define plot hieght based on the numer of pathways to be plotted
                             height=unit((selected_genes %in% heatmap.input$external_gene_name %>% sum())/2, "cm"),
                             width=unit(8, "cm"))
  
  heatmap
}

# remove intermediate files -----------------
rm(gene2name, kegg_t2g, pathway2gene, pathway2name, msigdb, biomart_species, kegg_species, species)










# ------------------------------

## Plot GSEA table
plotGSEAcompact <- function(.GSEA, title = "", cutoff = 0.05, subgroup = "all", truncate_label_at = 45) {
  .GSEA <- .GSEA %>% 
    filter(qvalues <= cutoff) %>%
    # Truncate long strings
    mutate(Description = str_trunc(Description, width = truncate_label_at, side = "right", ellipsis = ".."))
  
  # use it to set a common range scale
  max_NES <- .GSEA %>% filter(qvalues <= cutoff) %>% top_n(1, abs(NES)) %>% pull("NES") %>% abs()
  
  .GSEA.pos = .GSEA %>% filter(NES > 0)
  .GSEA.neg = .GSEA %>% filter(NES < 0)
  
  plot.pos <- .GSEA.pos %>%
    mutate(neglog10p = -log10(qvalues),
           Description = factor(Description, levels = 
                                  .GSEA.pos %>% arrange(NES) %>% pull(Description))) %>%
    ggplot(aes(x=Description, y=NES, fill=neglog10p, label=Description)) +
    geom_col() + 
    geom_text(y=0.1, hjust = 0, size = 3.5)+
    scale_y_continuous(limits = c(0, max_NES))+
    theme_bw() +
    scale_fill_continuous(type = "gradient",
                          high = "#DC3D2D", low = "#FED98B",
                          space = "Lab", na.value = "grey50",
                          label =  ~round(., digits = 1),
                          guide = "colourbar", name = "-Log10(adj p-value)")  +
    coord_flip() +
    ggtitle("GSEA: positive enrichment") +
    xlab("") +
    ylab("Normalized Enrichment Score") +
    theme(panel.grid = element_blank(),
          plot.title.position = "plot",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          title = element_text(size = 9))
  
  plot.neg <- .GSEA.neg %>%
    mutate(neglog10p = -log10(qvalues),
           Description = factor(Description, levels = 
                                  .GSEA.neg %>% arrange(desc(NES)) %>% pull(Description))) %>%
    ggplot(aes(x=Description, y=NES, fill=neglog10p, label=Description)) +
    geom_col() + 
    geom_text(y=-0.1, hjust = 1, size = 3.5)+
    scale_y_continuous(limits = c(-max_NES, 0))+
    theme_bw() +
    scale_fill_continuous(type = "gradient",
                          high = "#4A7AB7", low = "#C2E3EE",
                          space = "Lab", na.value = "grey50", guide = "colourbar",
                          label = ~round(., digits = 1),
                          name = "-Log10(adj p-value)")  +
    coord_flip() +
    ggtitle("GSEA: negative enrichment") +
    xlab("") +
    ylab("Normalized Enrichment Score") +
    theme(panel.grid = element_blank(),
          plot.title.position = "plot",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          title = element_text(size = 9)
          )
  
  # Fixed dimensions based on number of pathways
    plot.pos <-  set_panel_size(plot.pos,
                            # define plot hieght based on the numer of pathways to be plotted
                            height=unit((.GSEA.pos %>% filter(qvalues <= cutoff ) %>% nrow())/1.5, "cm"),
                            width=unit(12, "cm"))
    plot.neg <-  set_panel_size(plot.neg,
                                # define plot hieght based on the numer of pathways to be plotted
                                height=unit((.GSEA.neg %>% filter(qvalues <= cutoff ) %>% nrow())/1.5, "cm"),
                                width=unit(12, "cm"))
    
    
    if (subgroup == "pos") {
      .plot <- grid.arrange(plot.pos, ncol=1, top=grid::textGrob(title),
                            # Add 2 to free up some space for the title
                            heights =c((.GSEA.pos %>% filter(qvalues <= cutoff ) %>% nrow()/1.5) + 3))
    } else if (subgroup == "neg") {
      .plot <- grid.arrange(plot.neg, ncol=1, top=grid::textGrob(title),
                            # Add 2 to free up some space for the title
                            heights =c((.GSEA.neg %>% filter(qvalues <= cutoff ) %>% nrow()/1.5) + 3))
    } else if (subgroup == "all") {
      # plot the plot
      .plot <- grid.arrange(plot.pos, plot.neg, ncol=1, top=grid::textGrob(title),
                            # Add 2 to free up some space for the title
                            heights =c((.GSEA.pos %>% filter(qvalues <= cutoff ) %>% nrow()/1.5) + 3,
                                       .GSEA.neg %>% filter(qvalues <= cutoff ) %>% nrow()/1.5) + 2)
    }
    plot(.plot)
    # return a list, first item is the plot, second are the dimensions for saving it 
    list("plot" = .plot,
         "height" = case_when(
           # Add some margins for axis and others
           subgroup == "pos" ~ (.GSEA.pos %>% filter(qvalues <= cutoff) %>% nrow()/1.5 + 3),
           subgroup == "neg" ~ (.GSEA.neg %>% filter(qvalues <= cutoff) %>% nrow()/1.5 + 3),
           subgroup == "all" ~ (.GSEA %>% filter(qvalues <= cutoff) %>% nrow()/1.5 + 6)
         ))
}

