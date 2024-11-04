import pyBigWig
import os

# select chromosome
chr_selected = "chr14"

# Directories
input_dir = "/path/to/your/InputDir"  # Update this with the path to your input directory
output_dir = "/path/to/your/OutputDir"  # Update this with the path to your output directory

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Loop over each file in the input directory
for filename in os.listdir(input_dir):
    if filename.endswith(".bw"):
        # Construct full input file path
        input_bw = os.path.join(input_dir, filename)
        
        # Construct output file path
        output_bw = os.path.join(output_dir, filename.replace(".bw", "-" + chr_selected + ".bw"))

        # Open the input BigWig file
        bwIn = pyBigWig.open(input_bw)

        # Get the length of the selected chromosome
        chr_length = bwIn.chroms(chr_selected)

        # Open the output BigWig file for writing
        with pyBigWig.open(output_bw, 'w') as bwOutput:
            # Add header information (just for selected chromosome)
            bwOutput.addHeader([(chr_selected, chr_length)])

            # Write the data for the chromosome
            for x in bwIn.intervals(chr_selected, 0, chr_length):
                bwOutput.addEntries([chr_selected], [x[0]], ends=[x[1]], values=[x[2]])

        print(f"Processed {filename} -> {os.path.basename(output_bw)}")

        # Close the input BigWig file
        bwIn.close()

print("All files processed!")