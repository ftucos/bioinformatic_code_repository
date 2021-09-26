#!/bin/bash
set -ue -o pipefail
WD="/Volumes/TucosHDD/Bioinformatics/resources/genomes/homosapiens/STAR_RSEM/GENCODE_GRCh38_Index"

FASTA="/Volumes/TucosHDD/Bioinformatics/resources/genomes/homosapiens/GENCODE/GRCh38/GRCh38.primary_assembly.genome.fa"
GTF="/Volumes/TucosHDD/Bioinformatics/resources/genomes/homosapiens/GENCODE/GRCh38/gencode.v38.primary_assembly.annotation.gtf"
# need to provide the directory where STAR is located, not STAR executable itself
STAR="/Volumes/TucosHDD/Bioinformatics/resources/STAR-2.7.9a/bin/MacOSX_x86_64/"

cd $WD

rsem-prepare-reference --gtf $GTF \
	--star --star-path $STAR \
	--star-sjdboverhang 100 \
	--num-threads 8 \
	$FASTA \
	GRCh38