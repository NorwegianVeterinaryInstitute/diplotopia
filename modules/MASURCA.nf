// Asssembly with masurca
process MASURCA {
    //conda "/cluster/projects/nn9305k/src/miniconda/envs/"
    //conda (params.enable_conda ? 'bioconda::chewbbaca=3.3.1' : null)
    
    //container 'quay.io/biocontainers/masurca:4.1.0--pl5321hb5bd705_1'
    // some problem for some reason - try to use masurca directly in the container /usr/local/bin/masurca

    label 'process_high_memory_time'

    tag "$meta.id" 
    
    input:
    tuple val(meta), val(R1), val(R2), val(longRead), val(extraOption)

    output: 
    tuple val(meta), path("${meta.id}_out"), emit: masurca_assembly
    path "versions.yml", emit : versions


    script:
    // SECTION define script outputs
    def log_file = "${meta.id}_masurca.log"
    def out_dir="${meta.id}_out" 


    // !SECTION

    // SECTION preparing optional run with long reads and extra options
    def ploidy_value = meta.ploidy ?: '1'
    def hybrid_data_block = ''
    def linked_mates_param = '1'
    def jf_size = (meta.genome_size as long) * 20
    def flye_assembly_param = '0'


    if (meta.has_longRead) {
        hybrid_data_block = "NANOPORE=${longRead}"
        linked_mates_param = '0'
        flye_assembly_param = '1'
    } 
    // will not be used for now 
    def extra_config_param = extraOption ? extraOption : ''
    // !SECTION

    // SECTION GENERATE MASURCA config file (Gstring)
    def config_content = """
DATA
PE=pe ${params.pe} ${R1} ${R2}
${hybrid_data_block}
END

PARAMETERS
EXTEND_JUMP_READS=${params.extend_jump_reads}
GRAPH_KMER_SIZE = ${params.graph_kmer_size}
USE_LINKING_MATES = ${linked_mates_param}
USE_GRID=0
GRID_ENGINE=SGE
GRID_QUEUE=all.q
GRID_BATCH_SIZE=500000000
LHE_COVERAGE=${params.lhe_coverage}
LIMIT_JUMP_COVERAGE=${params.limit_jump_coverage}
CA_PARAMETERS = cgwErrorRate=${params.cgwErrorRate}
CLOSE_GAPS=${params.close_gaps}
NUM_THREADS = ${task.cpus}
JF_SIZE = ${jf_size}
SOAP_ASSEMBLY=${params.soap_assembly}
FLYE_ASSEMBLY=${flye_assembly_param}
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

    bash assemble.sh > ${log_file} 2>&1

    # Moving all in subdirectory for clarity
    mkdir ${out_dir}
    find . -mindepth 1 ! -name "${out_dir}" ! -type l -exec mv {} ${out_dir}/ \;
    """
    // !SECTION
}
