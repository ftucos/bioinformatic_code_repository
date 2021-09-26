# Download dbSNP 

Download VCF file from NCBI and add chr prefix to chromosomes

The following code requires samtools and cvbio, you can install both of them with conda

```shell
conda install cvbio samtools
```



```bash

WD=/Volumes/TucosHDD/Bioinformatics/resources/dbSNP
cd $WD

# download dbSNP from NCBI ----------------------
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/All_20180418.vcf.gz
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/All_20180418.vcf.gz.tbi

# add chr prefix to vcf ----------------------
wget https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCh38_ensembl2gencode.txt

## cvbio ends up convertig a bgzip to gzip
cvbio UpdateContigNames \
    -i All_20180418.vcf.gz \
    -o All_20180418.gencode_named.vcf.gz \
    -m GRCh38_ensembl2gencode.txt \
    --comment-chars '#' \
    --columns 0 \
    --skip-missing false

## convert from gzip back to bgzip
gunzip -c All_20180418.gencode_named.vcf.gz | bgzip -c > All_20180418.gencode_named.vcf.bgz

# rename it as a .gz file
rm -f All_20180418.gencode_named.vcf.gz
mv All_20180418.gencode_named.vcf.bgz All_20180418.gencode_named.vcf.gz

# generate tabindex ---------------------------
tabix -p vcf All_20180418.gencode_named.vcf.gz
```

