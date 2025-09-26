// https://mitra.stanford.edu/kundaje/marinovg/oak/programs/spades-cloudspades-paper/assembler/dipspades_manual.html
// Enventually add usage haplocontigs
process DIPSPADES {
    //conda "/cluster/projects/nn9305k/src/miniconda/envs/"
    //conda (params.enable_conda ? 'bioconda::chewbbaca=3.3.1' : null)
    
    container 'docker://pegi3s/spades'

    label 'process_high_memory'

    tag "$meta.id" 
    
    input:
    tuple val(meta), path(R1), path(R2)
    path(long_reads), optional : true
    val(options), optional : true

    output: 
    tuple val(meta), path("${meta.id}/*"), emit: masurca_assembly
    // Check path how to give path("CA/") - need assembly and need all results in output ? 
    path "versions.yml", emit : versions

    script:
    def short_reads_param = "-1 ${R1} -2 ${R2}"
    def long_reads_param = long_reads ? "-r ${long_reads}" : ""
    def advanced_options = "--hap haplocontigs.fasta"
    def log_file = "${meta.id}_dispades.log"

    """
    dipspades.py -v > versions.yml

    # To perform haplotype assembly of diploid genome - Illumina
    dipspades.py ${short_reads_param} \\
                ${advanced_options} \\
                -o ${meta.id}


    masurca -t $task.cpus \\
 
        ${long_reads_param} \\
        2>&1 | tee ${log_file} 
    """

}

