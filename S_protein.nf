Channel
    .fromPath('data/S_protein_samples_templates.csv')
    .splitCsv(header:true)
    .map{ row -> tuple(row.ID, file(row.sample), file(row.template))}
    .set{ mapping_ch }

process minimap2 {                                                              
    publishDir "./output/S_protein", mode: 'copy', overwrite: true                               
                                                                                
    input:                                                                      
    set ID, file(sample), file(template) from mapping_ch
                                                                                
    output:                                                                     
    file "*.sam" into sam_ch                                                    
                                                                                
    script:                                                                     
    """                                                                         
    minimap2 -t 8 -ax splice --sam-hit-only --cs=long \
      $template $sample > ${ID}.sam   
    """                                                                         
}                                                                               
                                                                                
process samtools_sort {                                                         
    publishDir "./output/S_protein", mode: 'copy', overwrite: true                               
                                                                                
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
    publishDir "./output/S_protein", mode: 'copy', overwrite: true                               
                                                                                
    input:                                                                      
    file bam from bam_ch                                                        
                                                                                
    output:                                                                     
    file "*.bam.bai"                                                            
                                                                                
    script:                                                                     
    """                                                                         
    samtools index $bam                                                         
    """                                                                         
}

