# RNAseq - Variant Calling

## 1. fastQC: Evaluate FASTQ quality

```bash
#!/bin/bash
# exit on error or unset variable, output all the command executed on the current shell
set -uex

WD=~/project/dir
FILE_DIR=$WD/data/fastq
OUT_DIR=$WD/logs/raw_fastqc
mkdir -p $OUT_DIR

find $FILE_DIR -name "*.fastq.gz" | parallel "fastqc {} -o $OUT_DIR"
```

## 2. Cutadapt: trim adapters

Trim adapters and evaluate new fastqc report

```bash
#!/bin/bash
# exit on error or unset variable, output all the command executed on the current shell
set -uex

WD=~/project/dir
FILE_DIR=$WD/data/fastq
OUT_DIR=$WD/processed/trimmed_fastq
OUT_FASTQC=$WD/logs/raw_fastqc
mkdir -p $OUT_DIR

# Copy files that don't require trimming ---------------------------------------------------
# NOTE: iterating on multiline string works only in bash and not in zsh shell
TO_LINK="
Sample_SW480_bulk/SW480_bulk_AGTTCC_L005_R1_001.fastq.gz
Sample_SW480_bulk/SW480_bulk_AGTTCC_L005_R2_001.fastq.gz
Sample_SW480_bulk1/SW480_bulk1_CGATGT_L004_R2_001.fastq.gz
Sample_SW480_bulk1/SW480_bulk1_CGATGT_L004_R1_001.fastq.gz
Sample_SW480_CD44highEpCAMlow/SW480_CD44highEpCAMlow_CTTGTA_L007_R2_001.fastq.gz
Sample_SW480_CD44highEpCAMlow/SW480_CD44highEpCAMlow_CTTGTA_L007_R1_001.fastq.gz
Sample_SW480_spheres/SW480_spheres_GCCAAT_L002_R1_001.fastq.gz
Sample_SW480_spheres/SW480_spheres_GCCAAT_L002_R2_001.fastq.gz"

for i in $TO_LINK
do
  cp $FILE_DIR/$i $OUT_DIR/$(basename -- $i)
done

# Trim read pairs ----------------------------------------------------------------------------------

# Adapter sequence obtained from  https://gist.github.com/photocyte/3edd9401d0b13476e60f8b104c2575f8
# Cut TruSeq Adapter: Index extracted from R1 file name & Illumina Single End PCR Primer 1 (reverse complemented) from R2

# TruSeq Adapter
ADAPTER_pt1="GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
ADAPTER_pt2="ATCTCGTATGCCGTCTTCTGCTTG"
#Illumina Single End PCR Primer 1 (reverse complemented) 
PRIMER="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"

# cutadapt can't work in multicore on python 3.8 on macOS so you need a conda enviroment with python 3.7 and cutadapt installed

R1=(
"Sample_SW480_bulk/SW480_bulk_ACAGTG_L007_R1_001.fastq.gz"
"Sample_SW480_bulk2/SW480_bulk2_CCGTCC_L004_R1_001.fastq.gz"
"Sample_SW480_CD44highEpCAMhigh/SW480_CD44highEpCAMhigh_CAGATC_L006_R1_001.fastq.gz"
"Sample_SW480_CD44highEpCAMhigh/SW480_CD44highEpCAMhigh_GCCAAT_L007_R1_001.fastq.gz"
"Sample_SW480_CD44highEpCAMhigh/SW480_CD44highEpCAMhigh_TGACCA_L004_R1_001.fastq.gz"
"Sample_SW480_CD44highEpCAMhigh2/SW480_CD44highEpCAMhigh2_GTCCGC_L004_R1_001.fastq.gz"
"Sample_SW480_CD44highEpCAMlow/SW480_CD44highEpCAMlow_ACAGTG_L002_R1_001.fastq.gz"
"Sample_SW480_CD44highEpCAMlow/SW480_CD44highEpCAMlow_CAGATC_L004_R1_001.fastq.gz"
"Sample_SW480_CD44highEpCAMlow/SW480_CD44highEpCAMlow_TGACCA_L006_R1_001.fastq.gz"
"Sample_SW480_spheres/SW480_spheres_AGTCAA_L007_R1_001.fastq.gz"
"Sample_SW480_spheres/SW480_spheres_AGTTCC_L004_R1_001.fastq.gz"
"Sample_SW480_spheres/SW480_spheres_CGATGT_L006_R1_001.fastq.gz"
)

R2=(
"Sample_SW480_bulk/SW480_bulk_ACAGTG_L007_R2_001.fastq.gz"
"Sample_SW480_bulk2/SW480_bulk2_CCGTCC_L004_R2_001.fastq.gz"
"Sample_SW480_CD44highEpCAMhigh/SW480_CD44highEpCAMhigh_CAGATC_L006_R2_001.fastq.gz"
"Sample_SW480_CD44highEpCAMhigh/SW480_CD44highEpCAMhigh_GCCAAT_L007_R2_001.fastq.gz"
"Sample_SW480_CD44highEpCAMhigh/SW480_CD44highEpCAMhigh_TGACCA_L004_R2_001.fastq.gz"
"Sample_SW480_CD44highEpCAMhigh2/SW480_CD44highEpCAMhigh2_GTCCGC_L004_R2_001.fastq.gz"
"Sample_SW480_CD44highEpCAMlow/SW480_CD44highEpCAMlow_ACAGTG_L002_R2_001.fastq.gz"
"Sample_SW480_CD44highEpCAMlow/SW480_CD44highEpCAMlow_CAGATC_L004_R2_001.fastq.gz"
"Sample_SW480_CD44highEpCAMlow/SW480_CD44highEpCAMlow_TGACCA_L006_R2_001.fastq.gz"
"Sample_SW480_spheres/SW480_spheres_AGTCAA_L007_R2_001.fastq.gz"
"Sample_SW480_spheres/SW480_spheres_AGTTCC_L004_R2_001.fastq.gz"
"Sample_SW480_spheres/SW480_spheres_CGATGT_L006_R2_001.fastq.gz"
)

for ((i = 0; i < ${#R1[@]} && i < ${#R2[@]}; i++))
do
  # extract indexes
  INDEX1=`sed 's/.*\([ACTG]\{6\}\).*/\1/' <<< ${R1[i]}`
  INDEX2=`sed 's/.*\([ACTG]\{6\}\).*/\1/' <<< ${R2[i]}`
  # verify indexes from R1 and R2 matches
  if [ "$INDEX1" != "$INDEX2" ]; then
    echo "Error: non matching Indexes between R1 and R2"
    exit 1
  fi
  cutadapt -a $ADAPTER_pt1$INDEX1$ADAPTER_pt2 -A $PRIMER \
    --minimum-length 20 \
    --cores 40 \
    --output $OUT_DIR/$(basename -- ${R1[i]}) \
    --paired-output $OUT_DIR/$(basename -- ${R2[i]}) \
    $FILE_DIR/${R1[i]} $FILE_DIR/${R2[i]} | tee -a $WD/logs/cutadapt.log
done

find $FILE_DIR -name "*.fastq.gz" | parallel "fastqc {} -o $OUT_FASTQC"
```

## 3. STAR: alignment

### Download Reference Genome

Be sure to have samtools and picard installed and that they are available on avery new terminal session automatically.

Otherwise you have to specify the path to the tools or update the PATH/load conda in the script `source ~/.bashrc`

```bash
#!/bin/bash
# exit on error or unset variable, output all the command executed on the current shell
set -uex

WD=/mnt/data/ccbc_environment/users/ftucci/resources/genomes/homosapiens
mkdir -p $WD/GENCODE_35/GRCh38

cd $WD/GENCODE_35/GRCh38
# Download Genome sequence (GRCh38.p13)
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/GRCh38.p13.genome.fa.gz
# Download Comprehensive gene annotation
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.chr_patch_hapl_scaff.annotation.gtf.gz

# Unzip
gunzip *.gz

#Generate Index
samtools faidx *.fa
# Generate Dict
picard CreateSequenceDictionary -R *.fa
```

### Build Genome Index

To improve speed at the alignment step sjdbOverhang should be set as readLength -1. In case of longer sequencing reads (more than >100) you should se it to the proper length

```bash
#!/bin/bash
# exit on error or unset variable, output all the command executed on the current shell
set -uex

STAR=/mnt/data/ccbc_environment/software/general/STAR-2.7.6a/STAR
REF_GENOME=~/userData/resources/genomes/homosapiens/GENCODE_35/GRCh38/GRCh38.p13.genome.fa
GENOME_ANNOTATIONS=~/userData/resources/genomes/homosapiens/GENCODE_35/GRCh38/GRCh38.p13.genome.fa/gencode.v35.chr_patch_hapl_scaff.annotation.gtf
OUT_STAR_INDEX=~/userData/resources/genomes/homosapiens/STAR/GRCh38

mkdir -p $OUT_STAR_INDEX
rm -rf $OUT_STAR_INDEX/*

cd $OUT_STAR_INDEX

$STAR --runMode genomeGenerate --runThreadN 96 --genomeDir $OUT_STAR_INDEX \
  --genomeFastaFiles $REF_GENOME \
  --sjdbGTFfile $GENOME_ANNOTATIONS \
  --sjdbOverhang 100
```

### Star Align

```bash
#!/bin/bash
# exit on error or unset variable, output all the command executed on the current shell
set -uex

STAR=/mnt/data/ccbc_environment/software/general/STAR-2.7.6a/STAR
GENOME_DIR=~/userData/resources/genomes/homosapiens/STAR/GRCh38


WD=~/userData/heterogeneity_SW480
FILE_DIR=$WD/processed/trimmed_fastq
OUTPUT_DIR=$WD/processed/STAR_alignment

mkdir -p $OUTPUT_DIR

# required only on MacOS
# increase the limit of open files from 256 to 50000 (required for sorting the BAM file)
# ulimit -n 5000

# Process Bulk ------------------------------------------------------------------------------------------------

# Parameter configuration as suggested in ENCODE example and Paper: Two-pass alignment improves novel splice junction quantification
# --outSAMmapqUnique set to 60 (default 255) to make quality scale compatible with GATK downstraam analysis

# Read group infos filled following https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups and https://support.sentieon.com/appnotes/read_groups/
# with unique @RG ID

$STAR --runThreadN 96 --genomeDir $GENOME_DIR \
  --readFilesIn $FILE_DIR/SW480_bulk_ACAGTG_L007_R1_001.fastq.gz,$FILE_DIR/SW480_bulk_AGTTCC_L005_R1_001.fastq.gz,$FILE_DIR/SW480_bulk1_CGATGT_L004_R1_001.fastq.gz,$FILE_DIR/SW480_bulk2_CCGTCC_L004_R1_001.fastq.gz \
  $FILE_DIR/SW480_bulk_ACAGTG_L007_R2_001.fastq.gz,$FILE_DIR/SW480_bulk_AGTTCC_L005_R2_001.fastq.gz,$FILE_DIR/SW480_bulk1_CGATGT_L004_R2_001.fastq.gz,$FILE_DIR/SW480_bulk2_CCGTCC_L004_R2_001.fastq.gz \
  --readFilesCommand gunzip -c \
  --outSAMtype BAM SortedByCoordinate \
  --twopassMode Basic \
  --outSAMattrRGline ID:SW480-BULK.C8NE5.7.ACAGTG PU:C8NE5ACXX.7 SM:SW480-BULK LB:BULK-1 PL:ILLUMINA , ID:SW480-BULK.C8NE5.5.AGTTCC PU:C8NE5ACXX.5 SM:SW480-BULK LB:BULK-2 PL:ILLUMINA , ID:SW480-BULK.C8NE5.4.CGATGT PU:C8NE5ACXX.4 SM:SW480-BULK LB:BULK-3 PL:ILLUMINA , ID:SW480-BULK.C8NE5.4.CCGTCC PU:C8NE5ACXX.4 SM:SW480-BULK LB:BULK-4 PL:ILLUMINA \
  --outSAMmapqUnique 60 \
  --alignSJDBoverhangMin 1 --scoreGenomicLengthLog2scale 0 --alignEndsType EndToEnd \
  --outFileNamePrefix $OUTPUT_DIR/SW480_bulk_

# Process CD44 Hi EpCAM Hi ------------------------------------------------------------------------------------
$STAR --runThreadN 96 --genomeDir $GENOME_DIR \
  --readFilesIn $FILE_DIR/SW480_CD44highEpCAMhigh_CAGATC_L006_R1_001.fastq.gz,$FILE_DIR/SW480_CD44highEpCAMhigh_GCCAAT_L007_R1_001.fastq.gz,$FILE_DIR/SW480_CD44highEpCAMhigh_TGACCA_L004_R1_001.fastq.gz,$FILE_DIR/SW480_CD44highEpCAMhigh2_GTCCGC_L004_R1_001.fastq.gz \
  $FILE_DIR/SW480_CD44highEpCAMhigh_CAGATC_L006_R2_001.fastq.gz,$FILE_DIR/SW480_CD44highEpCAMhigh_GCCAAT_L007_R2_001.fastq.gz,$FILE_DIR/SW480_CD44highEpCAMhigh_TGACCA_L004_R2_001.fastq.gz,$FILE_DIR/SW480_CD44highEpCAMhigh2_GTCCGC_L004_R2_001.fastq.gz \
  --readFilesCommand gunzip -c \
  --outSAMtype BAM SortedByCoordinate \
  --twopassMode Basic \
  --outSAMattrRGline ID:SW480-EPCAM_HI.C8NE5.6.CAGATC PU:C8NE5ACXX.6 SM:SW480-EPCAM_HI LB:EPCAM_HI-1 PL:ILLUMINA , ID:SW480-EPCAM_HI.C8NE5.7.GCCAAT PU:C8NE5ACXX.7 SM:SW480-EPCAM_HI LB:EPCAM_HI-2 PL:ILLUMINA , ID:SW480-EPCAM_HI.C8NE5.4.TGACCA PU:C8NE5ACXX.4 SM:SW480-EPCAM_HI LB:EPCAM_HI-3 PL:ILLUMINA , ID:SW480-EPCAM_HI.C8NE5.4.GTCCGC PU:C8NE5ACXX.4 SM:SW480-EPCAM_HI LB:EPCAM_HI-4 PL:ILLUMINA \
  --outSAMmapqUnique 60 \
  --alignSJDBoverhangMin 1 --scoreGenomicLengthLog2scale 0 --alignEndsType EndToEnd \
  --outFileNamePrefix $OUTPUT_DIR/SW480_CD44highEpCAMhigh_

# Process CD44 Hi EpCAM Low -----------------------------------------------------------------------------------
$STAR --runThreadN 96 --genomeDir $GENOME_DIR \
  --readFilesIn $FILE_DIR/SW480_CD44highEpCAMlow_ACAGTG_L002_R1_001.fastq.gz,$FILE_DIR/SW480_CD44highEpCAMlow_CAGATC_L004_R1_001.fastq.gz,$FILE_DIR/SW480_CD44highEpCAMlow_CTTGTA_L007_R1_001.fastq.gz,$FILE_DIR/SW480_CD44highEpCAMlow_TGACCA_L006_R1_001.fastq.gz \
  $FILE_DIR/SW480_CD44highEpCAMlow_ACAGTG_L002_R2_001.fastq.gz,$FILE_DIR/SW480_CD44highEpCAMlow_CAGATC_L004_R2_001.fastq.gz,$FILE_DIR/SW480_CD44highEpCAMlow_CTTGTA_L007_R2_001.fastq.gz,$FILE_DIR/SW480_CD44highEpCAMlow_TGACCA_L006_R2_001.fastq.gz \
  --readFilesCommand gunzip -c \
  --outSAMtype BAM SortedByCoordinate \
  --twopassMode Basic \
  --outSAMattrRGline ID:SW480-EPCAM_LOW.C8NE5.2.ACAGTG PU:C8NE5ACXX.2 SM:SW480-EPCAM_LOW LB:EPCAM_LOW-1 PL:ILLUMINA , ID:SW480-EPCAM_LOW.C8NE5.4.CAGATC PU:C8NE5ACXX.4 SM:SW480-EPCAM_LOW LB:EPCAM_LOW-2 PL:ILLUMINA , ID:SW480-EPCAM_LOW.C8NE5.7.CTTGTA PU:C8NE5ACXX.7 SM:SW480-EPCAM_LOW LB:EPCAM_LOW-3 PL:ILLUMINA , ID:SW480-EPCAM_LOW.C8NE5.6.TGACCA PU:C8NE5ACXX.6 SM:SW480-EPCAM_LOW LB:EPCAM_LOW-1 PL:ILLUMINA \
  --outSAMmapqUnique 60 \
  --alignSJDBoverhangMin 1 --scoreGenomicLengthLog2scale 0 --alignEndsType EndToEnd \
  --outFileNamePrefix $OUTPUT_DIR/SW480_CD44highEpCAMlow_

# Process Spheres ---------------------------------------------------------------------------------------------
$STAR --runThreadN 96 --genomeDir $GENOME_DIR \
  --readFilesIn $FILE_DIR/SW480_spheres_AGTCAA_L007_R1_001.fastq.gz,$FILE_DIR/SW480_spheres_AGTTCC_L004_R1_001.fastq.gz,$FILE_DIR/SW480_spheres_CGATGT_L006_R1_001.fastq.gz,$FILE_DIR/SW480_spheres_GCCAAT_L002_R1_001.fastq.gz \
  $FILE_DIR/SW480_spheres_AGTCAA_L007_R2_001.fastq.gz,$FILE_DIR/SW480_spheres_AGTTCC_L004_R2_001.fastq.gz,$FILE_DIR/SW480_spheres_CGATGT_L006_R2_001.fastq.gz,$FILE_DIR/SW480_spheres_GCCAAT_L002_R2_001.fastq.gz \
  --readFilesCommand gunzip -c \
  --outSAMtype BAM SortedByCoordinate \
  --twopassMode Basic \
  --outSAMattrRGline ID:SW480-SPHERES.C8NE5.7.AGTCAA PU:C8NE5ACXX.7 SM:SW480-SPHERES LB:SPHERES-1 PL:ILLUMINA , ID:SW480-SPHERES.C8NE5.4.AGTTCC PU:C8NE5ACXX.4 SM:SW480-SPHERES LB:SPHERES-1 PL:ILLUMINA , ID:SW480-SPHERES.C8NE5.6.CGATGT PU:C8NE5ACXX.6 SM:SW480-SPHERES LB:SPHERES-1 PL:ILLUMINA , ID:SW480-SPHERES.C8NE5.2.GCCAAT PU:C8NE5ACXX.2 SM:SW480-SPHERES LB:SPHERES-1 PL:ILLUMINA \
  --outSAMmapqUnique 60 \
  --alignSJDBoverhangMin 1 --scoreGenomicLengthLog2scale 0 --alignEndsType EndToEnd \
  --outFileNamePrefix $OUTPUT_DIR/SW480_spheres_

# Generate index file -----------------------------------------------------------------------------------------
cd $OUTPUT_DIR && parallel samtools index ::: *.bam

# Reorganize files in the proper directories ------------------------------------------------------------------
cd $OUTPUT_DIR

# Move Logs files
mkdir -p $WD/logs/STAR/
for i in *_Log*
do
  mv $i $WD/logs/STAR/$i
done

# Move BAM Files
mkdir -p $WD/processed/sandbox/bam
for i in *.out.bam{,.bai}
do
 mv $i $WD/processed/sandbox/bam/$i
done
```

## 4. BAM processing

```bash
#!/bin/bash
# exit on error or unset variable, output all the command executed on the current shell
set -uex

# Reference Files and Dirs
WD=~/userData/heterogeneity_SW480
FILE_DIR=$WD/processed/sandbox/bam/
REF_GENOME=~/userData/resources/genomes/homosapiens/GENCODE_35/GRCh38/GRCh38.p13.genome.fa
SNP1000G=/mnt/data/ccbc_environment/general/annotation/hg38/dbSNP/withChr_dbSNPv151_hg38_All_20180418.vcf.gz
GATK=/home/ftucci/userData/software/gatk-4.1.9.0/gatk

mkdir -p $WD/logs/picard
cd $FILE_DIR

# 1. Mark Duplicates ------------------------------------------------------------------------------------------
find ./ -name "*Aligned.sortedByCoord.out.bam" | sed "s/.*\///" | \
parallel picard MarkDuplicates -I ./{} -O $FILE_DIR{=s/_Aligned.sortedByCoord.out/.MarkDuplicates/=} --METRICS_FILE $WD/logs/picard/{=s/_Aligned.ByCoord.out.bam/MarkDuplicates.txt/=} -CREATE_INDEX true && \
# Remove STAR aligned raw files
rm -rf ./*Aligned.sortedByCoord.out.bam{,.bai}

# 2. Split reads spanning over junctions ----------------------------------------------------------------------
find ./ -name "*.MarkDuplicates.bam" | sed "s/.*\///" | \
parallel $GATK --java-options "-Xmx32G" SplitNCigarReads -R $REF_GENOME -I $FILE_DIR/{} -O $FILE_DIR/{=s/.MarkDuplicates.bam/.splitNtrim.bam/=} && \
# Remove temp files
rm -rf ./*.MarkDuplicates.bam{,.bai}

# 3. BQSR -----------------------------------------------------------------------------------------------------
find ./ -name "*.splitNtrim.bam" | sed "s/.*\///" | \
 parallel $GATK --java-options "-Xmx32G" BaseRecalibrator -R $REF_GENOME --known-sites $SNP1000G -I $FILE_DIR/{} -O $FILE_DIR/{=s/.splitNtrim.bam/.recal_data.table/=} && \
find ./ -name "*.splitNtrim.bam" | sed "s/.*\///" | \
  parallel $GATK --java-options "-Xmx32G" ApplyBQSR -R $REF_GENOME -I $FILE_DIR/{} --bqsr-recal-file $FILE_DIR/{=s/.splitNtrim.bam/.recal_data.table/=} -O $FILE_DIR/{=s/.splitNtrim.bam/.BQSR.bam/=}
rm -rf ./*.splitNtrim.bam{,.bai}

```

## ... Variant calling step



## InterVar

Note: in `--skip_annovar` mode, the out file should point to the annovar vcf file name without the ".hg38_multianno.vcf" appendix. The best input file to provide is the AVinput so InterVar don't have to generate that file by itself.

It's not possible to specify an output directory. The output file will be in the same directory of the source .txt file with the ".intervar" appendix

```bash
$INTERVAR/Intervar.py \
-t $INTERVAR/intervardb --input_type=AVinput --skip_annovar -b hg38 \
-i ./out/SW480_bulk_0.INDEL.avinput \
-o ./out/SW480_bulk_0.INDEL
```

**Expected output**

![image-20201105185519928](/Users/tucos/Library/Application Support/typora-user-images/image-20201105185519928.png)

