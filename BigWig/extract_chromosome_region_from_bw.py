import pyBigWig

# select chromosome
chr_selected = "chr14"

# Input BigWig file
input_bw = "/Volumes/TucosHDD/Bioinformatics/data/public_dataset/Reanalysis - Public BCa ChipSeq/Neyret-Kahn 2023 - Epigenomic mapping identifies an enhancer repertoire/data/Tumor Chip-seq/GSE193889/GSE193889/GSM5823194_HNNN129.bw"
# Output BigWig file
output_bw = "/Volumes/TucosHDD/Bioinformatics/data/public_dataset/Reanalysis - Public BCa ChipSeq/Neyret-Kahn 2023 - Epigenomic mapping identifies an enhancer repertoire/data/Tumor Chip-seq/GSE193889/GSE193889/GSM5823194_HNNN129-chr14.bw"

# Open the input BigWig file
bwIn = pyBigWig.open(input_bw)

# Get the length of the chromosome 
chr_length = bwIn.chroms()[chr_selected]

bwOutput = pyBigWig.open(output_bw,'w')
bwOutput.addHeader([(chr_selected,chr_length)]) # chromosome size

for x in bwIn.intervals(chr_selected,0,chr_length):
    bwOutput.addEntries([chr_selected],[x[0]],ends=[x[1]],values=[x[2]])

bwOutput.close()
bwIn.close()