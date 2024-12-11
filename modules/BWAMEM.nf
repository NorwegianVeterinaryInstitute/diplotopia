
// this is sam and not bam ! 
process BWAMEM {
	conda (params.enable_conda ? 'bioconda::bwa=0.7.8' : null)
	container 'quay.io/biocontainers/bwa:0.7.8--he4a0461_9'

    label 'process_high'
    
    tag "$ID" 

    input:
    tuple val(ID), path(assembly), path(R1), path(R2)

    output:
    tuple val(ID), path("*_temp.bam"), emit: shortbam_ch

    script:
    """
	bwa index $assembly
	bwa mem -t $task.cpus $assembly $R1 $R2 > ${ID}_temp.bam
    """
}

process BWAMEM_VAR {
	conda (params.enable_conda ? 'bioconda::bwa=0.7.8' : null)
	container 'quay.io/biocontainers/bwa:0.7.8--he4a0461_9'

    label 'process_high'
    
    tag "$ID" 

    input:
    tuple val(ID), path(R1), path(R2), path(assembly)

    output:
    tuple val(ID), path("*_map1.sam"), emit: shortbam_ch // will need to remove this one - in compass
    tuple val(ID), path("*_map1.sam"), emit: sam // check actually is not supposed to be compressed

    script:
    threads = task.cpus * 2
    """
	bwa index $assembly 
	bwa mem -t $threads -v 0 -R '@RG\tID:${ID}\tSM:${ID}' $params.bwamem_args  $assembly $R1 $R2 > ${ID}_map1.sam   
    """
}
// should I put all found alignments ? maybe make a pb with the variant calling
// removed the extra args -a -M

// TODO improve read group - now put all at sample ID ? `'@RG\tID:$ID\tSM:$SM\tLB:$LB\tPU:$PU\tPL:$PL'
//'@RG\tID:${ID}\tSM:${ID}' add PP:bwa (detected at samtools_merge_var) PP:bwa 
// https://stackoverflow.com/questions/72690191/bwa-mem-and-sambamba-read-group-line-error
// Tag PP:bwa not found in @PG records 
// [W::sam_hdr_link_pg] PG line with ID:samtools-5D89A8C2 has a PP link to missing program 'bwa'


// Need to find the way to write better 
process BWAMEM_VAR {
	conda (params.enable_conda ? 'bioconda::bwa=0.7.8' : null)
	container 'quay.io/biocontainers/bwa:0.7.8--he4a0461_9'

    label 'process_high'
    
    tag "$ID" 

    input:
    tuple val(ID), path(R1), path(R2), path(assembly)

    output:
    tuple val(ID), path("*_map1.sam"), emit: shortbam_ch // will need to remove this one - in compass
    tuple val(ID), path("*_map1.sam"), emit: sam // check actually is not supposed to be compressed

    script:
    threads = task.cpus * 2
    """
	bwa index $assembly 
	bwa mem -t $threads -v 0 -R '@RG\tID:${ID}\tSM:${ID}' $params.bwamem_args  $assembly $R1 $R2 > ${ID}_map1.sam   
    """
}