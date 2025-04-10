#!/bin/bash

# Read in a paths.txt file and a prefixes.txt file into arrays
# paths.txt is the directories of all pod5 files to be basecalled
# prefixes is like basename - the prefix to put on all output basecalled files
# THEY NEED TO BE IN THE SAME ORDER
mapfile -t paths < paths.txt
mapfile -t prefixes < prefixes.txt

# Check if arrays have same length
if [ ${#paths[@]} -ne ${#prefixes[@]} ]; then
    echo "Error: Number of paths does not equal prefixes"
    exit 1
fi

# Create output directory
mkdir -p basecalled

# Loop through both arrays simultaneously
for i in "${!paths[@]}"; do
    path="${paths[$i]}"
    prefix="${prefixes[$i]}"
    
    echo "Processing $prefix..."
    
		dorado basecaller \
		-v hac,inosine_m6A,m5C,pseU  \
		--min-qscore 10 \
		--emit-moves \
		--estimate-poly-a \
        --device "cuda:0,1" > "basecalled/${prefix}.bam"
        
    echo "Completed $prefix"
done
