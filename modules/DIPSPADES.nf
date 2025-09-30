// https://mitra.stanford.edu/kundaje/marinovg/oak/programs/spades-cloudspades-paper/assembler/dipspades_manual.html
// Enventually add usage haplocontigs
// you should set '--sc' flag if input data was obtained with MDA (single-cell) technology or --meta flag if processing metagenomic dataset)
// Memory limit (in Gb): 250 is per default

process DIPSPADES {
    //conda "/cluster/projects/nn9305k/src/miniconda/envs/"
    //conda (params.enable_conda ? 'bioconda::chewbbaca=3.3.1' : null)
    
    container 'docker://pegi3s/spades'

    label 'process_high_memory_time2'

    tag "$meta.id" 
    
    input:
    tuple val(meta), val(R1), val(R2), val(extraOption)

    output: 
    tuple val(meta),  path("${meta.id}_out"), emit: dipspades_assembly
    path "versions.yml", emit : versions

    script:
    def log_file = "${meta.id}_dispades.log"
    def out_dir="${meta.id}_out" 

    def short_reads_param = "-1 ${R1} -2 ${R2}"
    
    script:
    """
    dipspades.py -v > versions.yml

    # Haplotype assembly - short reads -diploid genome - Illumina
    dipspades.py ${short_reads_param} \\ 
                --memory ${task.memory} \\
                --threads ${task.cpus} \\
                ${extraOption} \\
                --cov-cutoff  auto \\
                -o ${meta.id} > ${log_file} 2>&1

    # Moving all in subdirectory for clarity
    mkdir ${out_dir}
    find . -mindepth 1 ! -name "${out_dir}" ! -type l -exec mv {} ${out_dir}/ \;
    """

}

