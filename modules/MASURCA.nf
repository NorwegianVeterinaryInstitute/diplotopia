// Asssembly with masurca
process MASURCA {
    //conda "/cluster/projects/nn9305k/src/miniconda/envs/"
    //conda (params.enable_conda ? 'bioconda::chewbbaca=3.3.1' : null)
    
    container 'quay.io/biocontainers/masurca:4.1.1--pl5321hb5bd705_0'

    label 'process_high_memory_time'

    tag "$meta.id" 
    
    input:
    tuple val(meta), path(R1), path(R2)
    val(longRead),  optional: true
    val(extraOption), optional: true
    val(ploidy_value)

    output: 
    tuple val(meta), path("*"), emit: masurca_assembly
    // will to rename path
    path "versions.yml", emit : versions

    script:
    def longRead_path = longRead ? file(longRead) : null 
    def log_file = "${meta.id}_masurca.log"
    def config_file_name = "${meta.id}_masurca.cfg"
    
    // Conditionally build the hybrid read line and flye parameter
    def hybrid_data_block = ''
    def flye_assembly_param = ''
    // will not be used for now 
    def extra_config_param = extraOption ? extraOption : ''


    // 'longRead' will be a non-null file path if provided, or null if not.
    // This Groovy check is correct.
    if (longRead) { 
        hybrid_data_block = "NANOPORE=${longRead}"
        flye_assembly_param = 'FLYE_ASSEMBLY=1'
    }

    // 2. Assemble the ENTIRE config file content into a single string (GString)
    def config_content = """
            DATA
            PE=pe 150 50 ${R1} ${R2}
            ${hybrid_data_block}
            END

            PARAMETERS
            USE_GRID=0
            EXTEND_JUMP_READS=0
            GRAPH_KMER_SIZE=auto
            LIMIT_JUMP_COVERAGE = 300
            LHE_COVERAGE=35
            MEGA_READS_ONE_PASS=0
            CLOSE_GAPS=1
            NUM_THREADS = \${task.cpus}
            JF_SIZE = 8240000000
            ${flye_assembly_param}
            ${extra_config_param} 
            END
    """

    // Execution 
    """
    masurca -v > versions.yml

    # Write the config to file
    echo """${config_content}""" > ${config_file_name}

    # Generate assemble.sh script
    masurca ${config_file_name}

    # Correct ploidy value in assemble.sh
    sed -i "s/echo 1 > PLOIDY.txt/echo ${ploidy_value} > PLOIDY.txt/" assemble.sh

    # Assembly
    bash assemble.sh
    """
}