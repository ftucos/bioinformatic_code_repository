# ANNOVAR: download reference databases

If you need the latest cosmic release, you need to build it by yourself. Follow Annovar Instructions

```bash
#!/bin/bash
ANNOVAR=/Volumes/TucosHDD/Bioinformatics/resources/annovar
cd $ANNOVAR

BUILD=hg38

# gene annotation from RefSeq 
$ANNOVAR/annotate_variation.pl -downdb -buildver $BUILD -webfrom annovar refGene humandb/ 

# gene annotation from Ensembl
$ANNOVAR/annotate_variation.pl -downdb -buildver $BUILD -webfrom annovar ensGene humandb/ 

# gene annotation from UCSC
$ANNOVAR/annotate_variation.pl -downdb -buildver $BUILD -webfrom annovar knownGene humandb/ 

# cytoBand annotation
$ANNOVAR/annotate_variation.pl --downdb --buildver $BUILD cytoBand humandb/

# alternative allele frequency (AAF) in the 1000 Genomes Project
$ANNOVAR/annotate_variation.pl --downdb --webfrom annovar --buildver $BUILD 1000g2015aug humandb/

# various functional deleteriousness prediction scores from dbNSFP (SIFT, CADD, MetaSVM etc.)
$ANNOVAR/annotate_variation.pl --downdb --webfrom annovar --buildver $BUILD dbnsfp35a humandb/

# variants reported in ClinVar (version 20200316)
$ANNOVAR/annotate_variation.pl --downdb --webfrom annovar --buildver $BUILD clinvar_20200316 humandb/

# ANNOVAR-compiled dbSNP (version 150)
$ANNOVAR/annotate_variation.pl --downdb --webfrom annovar --buildver $BUILD avsnp150 humandb/

# prediction of the splicing impact by Ada Boost and Random Forest
$ANNOVAR/annotate_variation.pl --downdb --webfrom annovar --buildver $BUILD dbscsnv11 humandb/

# minor allele frequenciy from gnomeAD project
$ANNOVAR/annotate_variation.pl --downdb --webfrom annovar --buildver $BUILD gnomad_genome humandb/

# alternative allele frequency in All subjects in the NHLBI-ESP project with 6500 exomes
$ANNOVAR/annotate_variation.pl --downdb --webfrom annovar --buildver $BUILD esp6500siv2_all humandb/

# Protein domain for variants
$ANNOVAR/annotate_variation.pl --downdb --webfrom annovar --buildver $BUILD dbnsfp31a_interpro humandb/

# RMSK (from UCSC)
curl http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/rmsk.txt.gz | gunzip -c > humandb/hg38_rmsk.txt
```

## Cosmic above v70

You have to download and build yourself the annotation file since Cosmic does not allows anymore the redistribution of it's files

```bash
# Cosmic Targeted Screen v92
wget http://www.openbioinformatics.org/annovar/download/prepare_annovar_user.pl
chmod +x prepare_annovar_user.pl

TEMP_DIR=$ANNOVAR/humandb/cosmic_temp
mkdir -p $TEMP_DIR

## key generated with key=$(echo "f.tucci@erasmusmc.nl:<your_password> | base64)
key="Zi50dWNjaUBlcmFzbXVzbWMubmw6U24tNzc3Nzk5OTkK"

## Send a download request for TSV file and extract the URL from the json response
TSV_URL=`curl -H "Authorization: Basic $key" https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v92/CosmicCompleteTargetedScreensMutantExport.tsv.gz | \
  sed | sed 's/{"url":"//' | sed 's/"}//'`

## Download and unzip the file
curl $TSV_URL | gunzip -c > $TEMP_DIR/CosmicCompleteTargetedScreensMutantExport.tsv

## Send a download request for VCF file and extract the URL from the json response
VCF_URL=`curl -H "Authorization: Basic $key" https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v92/VCF/CosmicCodingMuts.normal.vcf.gz | \
  sed | sed 's/{"url":"//' | sed 's/"}//'`

curl $VCF_URL | gunzip -c > $TEMP_DIR/CosmicCodingMuts.normal.vcf

$ANNOVAR/prepare_annovar_user.pl -dbtype cosmic $TEMP_DIR/CosmicCompleteTargetedScreensMutantExport.tsv -vcf $TEMP_DIR/CosmicCodingMuts.normal.vcf > humandb/hg38_cosmic92_targeted_coding.txt
```

