// Platanus 1.2.4 
// Platanus constructs consensus sequence of homologous chromosomes at first and tries to split into each haplotype sequence.
process PLATANUSS {
    //conda "/cluster/projects/nn9305k/src/miniconda/envs/"
    //conda (params.enable_conda ? 'bioconda::chewbbaca=3.3.1' : null)
    
    container 'cmonjeau/platanus:latest'
    

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

    Platanus assemble –o Pxut –f $GODOCKER_HOME/DRR02167[34]_[12].fastq –t 16 –m 128 2> assemble.log



    /root/src/redundans/redundans.py --version > versions.yml

    # important to put long_reads param at the end because create as string, in on new line produces error
    /root/src/redundans/redundans.py -v ${short_reads_param} ${long_reads_param} \\
                --mem ${task.mem} --tmp . \\
                -o ${meta.id} \\
                -t ${task.cpus} > ${log_file} 2>&1
    
    """
}



// PLACE Holder - we can try later 
// https://cell-innovation.nig.ac.jp/maser/ToolsList/P000001849.html
// Platanus-allee tries to construct each haplotype sequence from the beginning and pair them as homologous chromosomes, 
// process PLATANUS_ALLELE {

// }