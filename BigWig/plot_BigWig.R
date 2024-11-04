setwd("/Volumes/TucosHDD/Bioinformatics/workspace/BLADDER_cancer_project/bladder_bulkRNAseq_cell_lines_NUMB_hi_vs_low")
library(tidyverse)
library(patchwork)

NUMB_transcript <- c("ENST00000555238", "ENST00000356296", "ENST00000557597", "ENST00000554546")

# sample TPM ---------------------
ensembl2symbol <- biomaRt::getBM(mart = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl"),
                                 attributes = c("external_gene_name", "ensembl_gene_id", "ensembl_transcript_id", "gene_biotype", "transcript_biotype", "description")) %>%
  ## Remove source info in the description field
  mutate(description = str_remove(description, "\\s?\\[.*\\]\\s?"))

tpm <- read_tsv("/Volumes/TucosHDD/Bioinformatics/data/bladder_cancer/bulkRNAseq_cell_lines/star_rsem/rsem.merged.transcript_tpm.tsv") %>%
  rename("ensembl_gene_id" = "gene_id",
         "ensembl_transcript_id" = "transcript_id") %>%
  mutate(ensembl_gene_id = str_remove(ensembl_gene_id, "\\.[0-9]+$"),
         ensembl_transcript_id = str_remove(ensembl_transcript_id, "\\.[0-9]+$")) %>%
  gather(key="sample_id", value="TPM", -c(1:2)) %>%
  mutate(logTPM = log2(TPM+1))

st.err <- function(x) {sd(x, na.rm=TRUE)/sqrt(length(x[!is.na(x)]))}

logtpm.average <- tpm %>%
  mutate(sample_id = str_remove(sample_id, "_R[12]")) %>%
  group_by(ensembl_gene_id, ensembl_transcript_id, sample_id) %>%
  summarize(mean = mean(logTPM),
            st.err = st.err(logTPM))

logtpm.average <- logtpm.average %>%
  right_join(ensembl2symbol, .) %>%
  mutate(condition = case_when(
    sample_id == "5637" ~ "NUMB Low",
    sample_id == "CLS439" ~ "NUMB Low",
    sample_id == "HT1376" ~ "NUMB Low",
    sample_id == "KK47" ~ "NUMB High",
    sample_id == "RT112" ~ "NUMB High",
    sample_id == "RT4" ~ "NUMB High"),
    condition = factor(condition, levels = c("NUMB Low", "NUMB High")))

tpm %>%
  mutate(condition = case_when(
    str_detect(sample_id, "5637") ~ "NUMB Low",
    str_detect(sample_id, "CLS439") ~ "NUMB Low",
    str_detect(sample_id, "HT1376") ~ "NUMB Low",
    str_detect(sample_id, "KK47") ~ "NUMB High",
    str_detect(sample_id,"RT112") ~ "NUMB High",
    str_detect(sample_id, "RT4") ~ "NUMB High"),
    condition = factor(condition, levels = c("NUMB Low", "NUMB High"))) %>%
  left_join(ensembl2symbol) %>%
  filter(external_gene_name == "NUMB",
         transcript_biotype == "protein_coding",
         !ensembl_transcript_id %in% c("ENST00000554521", "ENST00000555307", "ENST00000555859", "ENST00000555987", "ENST00000557597", "ENST00000559312", "ENST00000560335", "ENST00000356296", "ENST00000555394", "ENST00000554818", "ENST00000555738", "ENST00000556772")) %>%
  #filter(ensembl_transcript_id %in% c("ENST00000555238", "ENST00000356296", "ENST00000557597", "ENST00000554546")) %>%
  mutate(
  #   transcript_id = case_when(
  #   ensembl_transcript_id == "ENST00000555238" ~ "NUMB1",
  #   ensembl_transcript_id == "ENST00000356296" ~ "NUMB2",
  #   ensembl_transcript_id == "ENST00000557597" ~ "NUMB3",
  #   ensembl_transcript_id == "ENST00000554546" ~ "NUMB4"
  # ),
  sample = str_remove(sample_id, "_R[12]")) %>%
  # keep the genes ordered by p-value and log2FC
  ggplot(aes(x=ensembl_transcript_id, y=logTPM, fill=condition, group=sample_id)) +
  geom_col(position = "dodge")+
  facet_wrap(condition~sample, scales="free_x") +
  xlab("") + ylab("Log(TPM+1)")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        legend.position = "none")

# Plot exon usage --------------------------------
# load required data
library(GenomicFeatures)
library(wiggleplotr)
library(tidyverse)

txdb <- makeTxDbFromGFF("/Volumes/TucosHDD/Bioinformatics/resources/genomes/homosapiens/GENCODE/GRCh38/gencode.v38.primary_assembly.annotation.gtf",
                        format="gtf" , organism='Homo sapiens', dataSource="GRCh38 GENCODE v38 - primary_assembly")

exons = exonsBy(txdb, by = "tx", use.names = TRUE)
names(exons) = str_remove(names(exons), "\\.[0-9]+$")
cdss = cdsBy(txdb, by = "tx", use.names = TRUE)
names(cdss) = str_remove(names(cdss), "\\.[0-9]+$")
annotations = biomaRt::getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "external_gene_name", "strand", "gene_biotype", "transcript_biotype"),
                             mart = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl"))%>%
  dplyr::rename("gene_id" = "ensembl_gene_id", 
         "transcript_id" = "ensembl_transcript_id",
         "gene_name" = "external_gene_name")

# ATAC input data
sample_data.ATAC = data.frame(
  sample_id = c("5637", "CLS439", "HT1376", "KK47", "RT112", "RT4"), 
  condition = factor(c("NUMB Low", "NUMB Low", "NUMB Low", "NUMB High", "NUMB High", "NUMB High"), levels = c("NUMB Low", "NUMB High")),
  # bigwig are prescaled
  scaling_factor = 1,
  bigWig = c(
    "processed/bam/mergedReplicate/bigwig/Sample_5637.NFR.bw",
    "processed/bam/mergedReplicate/bigwig/Sample_CLS439.NFR.bw",
    "processed/bam/mergedReplicate/bigwig/Sample_HT1376.NFR.bw",
    "processed/bam/mergedReplicate/bigwig/Sample_KK47.NFR.bw",
    "processed/bam/mergedReplicate/bigwig/Sample_RT112.NFR.bw",
    "processed/bam/mergedReplicate/bigwig/Sample_RT4.NFR.bw"
  )) %>%
  mutate(track_id = sample_id, colour_group = condition)

# scale RNAseq data -----------------------------
extract_scaling <- function(x) {data.frame(sample_id = str_extract(x, "(?<=deseq2_qc\\/size_factors\\/)[^\\.]+"),
                                           scaling_factor = read_tsv(x,col_names = F)$X1)}

scaling.factors <- list.files("/Volumes/TucosHDD/Bioinformatics/data/bladder_cancer/bulkRNAseq_cell_lines/star_rsem/deseq2_qc/size_factors",
                              pattern = ".txt", full.names = T) %>%
  map(extract_scaling) %>%
  do.call(what=rbind) %>%
  mutate(sample_id = str_replace(sample_id, pattern = "X5637_", replacement = "5637_"))

# exponentiate scaling since the size factors by deseq are computed on log2 values
scaling.factors.exponentiate <- scaling.factors %>%
  mutate(scaling_factor = 2^scaling_factor/2)

# RNAseq input data
sample_data.RNA = data.frame(
  bigWig = list.files("/Volumes/TucosHDD/Bioinformatics/data/bladder_cancer/bulkRNAseq_cell_lines/star_rsem/bigwig",
                      pattern =".reverse.bigWig", full.names = T)
) %>%
  mutate(sample = str_extract(bigWig, "(?<=star_rsem\\/bigwig\\/)[^_]+"),
         sample_id = str_extract(bigWig, "(?<=star_rsem\\/bigwig\\/)[^\\.]+"),
         condition = case_when(
           sample == "5637" ~ "NUMB Low",
           sample == "CLS439" ~ "NUMB Low",
           sample == "HT1376" ~ "NUMB Low",
           sample == "KK47" ~ "NUMB High",
           sample == "RT112" ~ "NUMB High",
           sample == "RT4" ~ "NUMB High"),
         condition = factor(condition, levels = c("NUMB Low", "NUMB High")),
         # scaling_factor = 1,
         track_id = sample_id,
         colour_group = condition) %>% left_join(scaling.factors.exponentiate)

# Plot
plot <- plotCoverage(exons = exons[NUMB_transcript], cdss = cdss[NUMB_transcript], transcript_annotations = annotations,
                     track_data = sample_data.RNA, 
                     coverage_type = 'both',
                     flanking_length = c(1000, 1000),
                     fill_palette = c("#0073C2", "#CD534C"), 
                     rescale_introns = T, transcript_label = F, return_subplots_list = T
)
((plot$coverage_plot +
    geom_vline(aes(xintercept = 2270), size = 0.5, color="#2D3E50") + 
    geom_vline(aes(xintercept = 2455), size = 0.5, color="#2D3E50") + 
    geom_vline(aes(xintercept = 3570), size = 0.5, color="#2D3E50") + 
    geom_vline(aes(xintercept = 3650), size = 0.5, color="#2D3E50") + 
  facet_grid(track_id~., scales="free_y")) /
(plot$tx_structure + 
    geom_vline(aes(xintercept = 2270), size = 0.5, color="#2D3E50") + 
    geom_vline(aes(xintercept = 2455), size = 0.5, color="#2D3E50") +
   geom_vline(aes(xintercept = 3570), size = 0.5, color="#2D3E50") + 
   geom_vline(aes(xintercept = 3650), size = 0.5, color="#2D3E50"))) +
  plot_layout(heights = c(5, 1)) +
  plot_annotation(title = "NUMB exon usage")

ggsave("output/NUMB_exon_usage.png", plot=last_plot(), height = 3000, width = 2500, units = "px", scale = 1.5)
  