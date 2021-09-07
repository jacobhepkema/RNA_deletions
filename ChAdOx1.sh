#!/bin/bash

echo "Aligning to plasmid"
minimap2 -t 8 -ax splice --sam-hit-only --cs=long templates/ChAdOx1_template.fa data/SRR13320597.fastq > output/ChAdOx1/SRR13320597_template.sam
echo "Converting to BAM"
samtools sort output/ChAdOx1/SRR13320597_template.sam > output/ChAdOx1/SRR13320597_template.bam
rm output/ChAdOx1/SRR13320597_template.sam
echo "Indexing BAM"
samtools index output/ChAdOx1/SRR13320597_template.bam
echo "Done"

# Add output script
