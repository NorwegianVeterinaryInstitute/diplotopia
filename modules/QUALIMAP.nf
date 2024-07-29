// Statistics from mapping
//path(*_coverage)
process QUALIMAP {
    conda "/cluster/projects/nn9305k/src/miniconda/envs/qualimap2"
    container 'quay.io/biocontainers/qualimap:2.2.2d--hdfd78af_2'

    label 'process_high'

    tag "$ID"

    input:
    tuple val(ID), path(bamfile)
    
    output: 
    path("${ID}"), emit: qualimap_ch

    script:
    // dynamically assign memory process - but not all memody
    javamem = "${task.memory.toGiga()-4}G"

    """
    qualimap bamqc \\
    --java-mem-size=${javamem} -nt $task.cpus \\
    -bam $bamfile --outdir . \\
    --paint-chromosome-limits \\
    --collect-overlap-pairs \\
    --output-genome-coverage _coverage
    
    mkdir ${ID}
    mv *_stats ${ID}
    mv *_coverage ${ID}
    """
}
// This creates an error and make pipeline bug find solution
// qualimap bamqc > qualimap.version