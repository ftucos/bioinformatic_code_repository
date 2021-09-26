# Download and setup reference genome

Download Reference genome GRCh38 for bioinformatic analysis from GENCODE https://www.gencodegenes.org/human/

GENCODE annotation is essentially identical to the Ensembl annotation (including way more transcript than the ones annotated by RefSeq). The main difference from Ensembl is that chromosomes includes the "chr" prefix.

Also,  genes which are common to the human chromosome X and Y PAR regions can be found twice in the GENCODE GTF, while they are shown only for chromosome X in the Ensembl file.

More info at: https://www.gencodegenes.org/pages/faq.html#:~:text=What%20is%20the%20difference%20between%20GENCODE%20GTF%20and%20Ensembl%20GTF,X%20in%20the%20Ensembl%20file.

```bash
#!/bin/bash

WD=/mnt/data/ccbc_environment/users/ftucci/resources/genomes/homosapiens
mkdir -p $WD/GENCODE_35/GRCh38

cd $WD/GENCODE_35/GRCh38
# Download Gencode GRCh38 primary genome assembly (chromosomes and scaffolds)
# most aligners don't deal with patches and haplotypes with the final result of simply increasing the number of multimappers
mkdir -p $GENOME_DIR
cd $GENOME_DIR
# Download fasta sequence
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/GRCh38.primary_assembly.genome.fa.gz
# Download annotation
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.primary_assembly.annotation.gtf.gz

# Gunzip
gunzip $GENOME_DIR/*.gz


# be sure you have Picard and samtools installed in the PATH

# Generate Index
samtools faidx $GENOME_DIR/*.fa
# Generate Dict
picard CreateSequenceDictionary -R $GENOME_DIR/*.fa
```



