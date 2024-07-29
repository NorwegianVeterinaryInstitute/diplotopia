// modules to work with sam and bam files 
// https://lomereiter.github.io/sambamba/docs/sambamba-view.html

// marking duplicates (or removing?) for variant calling with freebayes

// maybe add 
// tag "$ID" 

process SAMBAMBA_MARKDUP {
    //conda ""
    container 'quay.io/biocontainers/sambamba:1.0.1--h6f6fda4_1'

    label 'process_high'

    input:
    tuple val(ID), val(ploidy), path(bam)

    output:
    tuple val(ID), val(ploidy), path("*_markdup.bam"), emit: bam
    path("*")

    script:
    threads = task.cpus * 2
    if (params.markdup == "mark") {
        """
        sambamba markdup -t $threads -p $bam ${ID}_markdup.bam  2>&1 | tee ${ID}_sambamba_markup.log
        #sambamba -h > sambamba.version
        """
    } else if (params.markdup == "remove") {
        """
        sambamba markdup -t $threads -p -r $bam ${ID}_markdup.bam 2>&1 | tee ${ID}_sambamba_markup.log
        #sambamba -h > sambamba.version
        """
    } 
    
}

// exporting samba version makes it bug - find a solution 
// can be without label