#!/bin/bash

set -uex -o pipefail

WD="/Volumes/TucosHDD/Bioinformatics/data/Fodde_Lab/WES_reprocessing-OSCC/"
FILEDIR="$WD/result/variant_calling/controlfreec/"
OUTDIR="$WD/custom_results/controlfreec/"
REFBED="/Volumes/TucosHDD/Bioinformatics/resources/genomes/homosapiens/GENCODE/GRCh38/gencode.v38.primary_assembly.annotation.bed"

files=(
"$FILEDIR/CA1/CA1.p.value.txt"
"$FILEDIR/LM/LM.p.value.txt"
"$FILEDIR/LUC4/LUC4.p.value.txt"
"$FILEDIR/SCC4/SCC4.p.value.txt"
"$FILEDIR/SCC9/SCC9.p.value.txt"
"$FILEDIR/SCC15/SCC15.p.value.txt"
"$FILEDIR/SCC25/SCC25.p.value.txt"
)

mkdir -p $OUTDIR

# .p.value.txt is formattes as .bed as follows (header removed):
# chr	start	end	copy_number|status|genotype|uncertainty|WilcoxonRankSumTestPvalue|KolmogorovSmirnovPvalue


for file in ${files[@]}; do
	sample_name=$(echo $file | sed 's/\.p\.value\.txt//g' | sed 's/.*\///g')
	echo $sample_name
	awk -F '\t' 'NR > 1 { print "chr"$1, $2, $3, $4 "|" $5 "|" $6 "|" $7 "|" $8 "|" $9 }' $file | bedmap --echo --echo-map-id --delim '\t' --multidelim  '|' - $REFBED | sed -e 's/ /\t/g' > ${OUTDIR}/${sample_name}.annotated.p.value.txt
done
