# nf-core RNAseq processing

```bash
#!/usr/bin/bash
set -uex -o pipefail

WD=~/userData/workspace/siNUMB_RNAseq_processing/

# this pipeline requires the developmental version of nextflow
NEXTFLOW=~/userData/resources/tools/nextflow/nextflow-21.02.0-edge-all
GENOME_SEQUENCE=~/userData/resources/genomes/gencode/hsapiens/GRCh38/GRCh38.primary_assembly.genome.fa
GENOME_ANNOTATION=~/userData/resources/genomes/gencode/hsapiens/GRCh38/gencode.v36.primary_assembly.annotation.gtf
STAR_INDEX=~/userData/resources/genomes/STAR/GRCh38_overhang_50bp
RSEM_INDEX=~/userData/resources/genomes/RSEM/GRCh38_gencode

# downloa data from the bucket: gsutil -m cp -r gs://bladder_rna-seq/* .

cd $WD
# when using "star_rsem" aligner option, only the RSEM index field is parsed and you should include the STAR index withing the RSEM folder
$NEXTFLOW run nf-core/rnaseq \
	--input data/Bladder_siNUMB_samplesheet.csv \
	--fasta $GENOME_SEQUENCE \
	--gtf $GENOME_ANNOTATION --gencode \
	--rsem_index $RSEM_INDEX \
	--aligner 'star_rsem' \
	--skip_markduplicates \
	-profile docker \
	-resume backstabbing_bernard \
	--max_cpus 16 \
	--max_memory '100.GB' \
	--email_on_fail 'fran.tucci@gmail.com'
```

