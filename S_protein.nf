samples_file = file("data/S_protein_samples.txt")
samples_names = samples_file.readLines()
samples_ch = Channel.fromList(samples_names)
genome = file("data/hg38.fa")

process minimap2 {
    publishDir ".", mode: 'copy', overwrite: true

    input:
    path c_sample from samples_ch

    output:
    file "*.sam" into sam_ch

    script:
    """
    minimap2 -t 10 -ax splice -uf --sam-hit-only --cs=long \
      $genome $c_sample > ${c_sample.baseName}_hg38.sam
    """
}

process samtools_sort {
    publishDir ".", mode: 'copy', overwrite: true

    input:
    file sam from sam_ch
    
    output:
    file "*.bam" into bam_ch

    script:
    """
    samtools sort $sam > ${sam.baseName}.bam
    """
}

process samtools_index {
    publishDir ".", mode: 'copy', overwrite: true

    input:
    file bam from bam_ch

    output:
    file "*.bam.bai"

    script:
    """
    samtools index $bam
    """
}
