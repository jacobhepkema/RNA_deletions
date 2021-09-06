#!/bin/bash

echo "Aligning to plasmid"
minimap2 -t 8 -ax splice --sam-hit-only --cs=long templates/ALKBH5_mRuby3_plasmid.fa data/SRR9646144.fastq > output/ALKBH5/SRR9646144_plasmid.sam
echo "Converting to BAM"
samtools sort output/ALKBH5/SRR9646144_plasmid.sam > output/ALKBH5/SRR9646144_plasmid.bam
rm output/ALKBH5/SRR9646144_plasmid.sam
echo "Indexing BAM"
samtools index output/ALKBH5/SRR9646144_plasmid.bam
echo "Done"

# Add output script
