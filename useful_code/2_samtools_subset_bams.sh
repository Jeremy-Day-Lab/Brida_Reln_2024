#!/bin/bash

# Assuming you move/copy all BAMs from CellRanger output into working dir

# Subset to create BAMs with only reads from consolidated "Virus" chromosome (PX02070015)
mkdir virus_bams
i=0

for b in *.bam
do
    i=$((i+1))
    outfile="virus_bams/subset_virus_$i.bam"
    samtools view -h $b PX02070015 | samtools view -bS > $outfile
done

# Create BAM index files
cd virus_bams

for b in *.bam
do
    echo "Indexing: "$b
    samtools index $b
done

# Merge subset BAMs into one and index
samtools merge virus_bams.bam *.bam
samtools index virus_bams.bam

# Subset bams by cluster/cellType
mkdir cluster_bams
declare -a cluster_array=("astro" "d1" "d2_1" "d2_2" "d3" "gaba" "grm8" "micro" "mural" "olig" "poly" "pvalb" "sst")
echo ${cluster_array[@]}

# Iterate over clusters
for c in ${cluster_array[@]}
do
    # Iterate over GEM wells
    for x in 1 2 3 4 5 6 7 8
    do
        barcodefile="barcodes/${c}_${x}_barcodes.csv"
        # If barcode file exists for this cluster_GEM
        if [ -f $barcodefile ]
        then
            echo "Running subset for ${c}_${x}"
            outfile="cluster_bams/subset_virus_${c}_${x}.bam"
            # Change path to where you saved subset-bam_macos
            ~/subset-bam_macos --bam virus_bams.bam --cell-barcodes $barcodefile --out-bam $outfile
        # Skips iteration if barcode file doesn't exist for that cluster_GEM
        else
            echo "Skipping ${c}_${x}"
        fi
    done
done


