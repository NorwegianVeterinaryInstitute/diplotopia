// Asssembly with masurca
// https://github.com/aakashsur/docker-masurca/blob/main/Dockerfile
// lack ps but otherwise might be ok 
// Difficulties to make work in container 
process MASURCA {
    conda (params.enable_conda ? 'bioconda::masurca=4.1.1' : null)    
    // container 'quay.io/biocontainers/masurca:4.1.3--h6b3f7d6_1'
    // container 'docker://aakashsur/masurca'
    // 4.0.3 


    //stageInMode 'copy' 

    label 'process_high_memory_time'

    tag "$meta.id" 
    
    input:
    tuple val(meta), 
        path(R1, stageAs: 'R1.fastq.gz'), 
        path(R2, stageAs: 'R2.fastq.gz'), 
        path(longRead, stageAs: 'longRead.fa'), 
        val(extraOption)

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
        hybrid_data_block = "NANOPORE=longRead.fa"
        linked_mates_param = '0'
        flye_assembly_param = '1'
    } 
    // will not be used for now 
    def extra_config_param = extraOption ? extraOption : ''
    // !SECTION

    // SECTION GENERATE MASURCA config file (Gstring)
    def config_content = """
DATA
PE=pe ${params.pe} R1.fastq.gz R2.fastq.gz
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
    # FIXE need to find a proper container with masurca installed
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
