// Asssembly with masurca
process MASURCA {
    //conda "/cluster/projects/nn9305k/src/miniconda/envs/"
    //conda (params.enable_conda ? 'bioconda::chewbbaca=3.3.1' : null)
    
    //container 'quay.io/biocontainers/masurca:4.1.0--pl5321hb5bd705_1'

    label 'process_high_memory_time'

    tag "$meta.id" 
    
    input:
    tuple val(meta), val(R1), val(R2), val(longRead), val(extraOption), val(ploidy_value)

    output: 
    tuple val(meta), path("*"), emit: masurca_assembly
    path "versions.yml", emit : versions


    script:
    // SECTION define script outputs
    def log_file = "${meta.id}_masurca.log"

    // !SECTION

    // SECTION preparing optional run with long reads and extra options
    def hybrid_data_block = ''
    def flye_assembly_param = ''

    if (meta.has_longRead) {
        hybrid_data_block = "NANOPORE=${longRead}"
        flye_assembly_param = 'FLYE_ASSEMBLY=1'
    } 
    // will not be used for now 
    def extra_config_param = extraOption ? extraOption : ''
    // !SECTION

    // SECTION GENERATE MASURCA config file (Gstring)
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
NUM_THREADS = ${task.cpus}
JF_SIZE = 8240000000
${flye_assembly_param}
${extra_config_param} 
END
   """
   // !SECTION
  
  


    // SECTION MASURCA execution 
    """
    # FIXME could not make it work with containers I tried. Find way to use container. Missing gz ? 
    # maybe need to stage file but then pb when missing files ... and paths
    module load MaSuRCA/4.1.0-GCC-11.3.0
    masurca -v > versions.yml
    echo "${config_content}" > ${meta.id}_masurca.config

    masurca ${meta.id}_masurca.config

    sed -i "s/echo 1 > PLOIDY.txt/echo ${ploidy_value} > PLOIDY.txt/" assemble.sh

    bash assemble.sh

    """
}
