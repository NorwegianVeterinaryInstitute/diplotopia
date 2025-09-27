// https://github.com/Gabaldonlab/redundans/tree/master#
// Contains platanus assembler 


// Asssembly with masurca
process REDUNDANS {
    //conda "/cluster/projects/nn9305k/src/miniconda/envs/"
    //conda (params.enable_conda ? 'bioconda::chewbbaca=3.3.1' : null)
    
    container 'cgenomics/redundans:latest'

    label 'process_high_memory_time'

    tag "$meta.id" 
    
    input:
    tuple val(meta), val(R1), val(R2), val(longRead)

    output: 
    tuple val(meta), path("*"), emit: redundans_assembly
    path "versions.yml", emit : versions


    script:
    // SECTION define script outputs
    def log_file = "${meta.id}_redundans.log"

    def short_reads_param = "-i ${R1} ${R2}"
    def long_reads_param = '' 
    if (meta.has_longRead) {
        long_reads_param = "--longreads ${longRead}"
    }

    // SECTION MASURCA execution 
    """
    # FIXME could not make it work with containers I tried. Find way to use container. Missing gz ? 
    # maybe need to stage file but then pb when missing files ... and paths

    /root/src/redundans/redundans.py --version > versions.yml

    /root/src/redundans/redundans.py -v ${short_reads_param} ${long_reads_param} \\
                -o ${meta.id} \\
                -t ${task.cpus} 2>&1 | tee ${log_file}
    
    """
}
