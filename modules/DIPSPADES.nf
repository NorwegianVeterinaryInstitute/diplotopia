// https://mitra.stanford.edu/kundaje/marinovg/oak/programs/spades-cloudspades-paper/assembler/dipspades_manual.html
// Enventually add usage haplocontigs
process DIPSPADES {
    //conda "/cluster/projects/nn9305k/src/miniconda/envs/"
    //conda (params.enable_conda ? 'bioconda::chewbbaca=3.3.1' : null)
    
    container 'docker://pegi3s/spades'

    label 'process_high_memory'

    tag "$meta.id" 
    
    input:
    tuple val(meta), val(R1), val(R2), val(extraOption)

    output: 
    tuple val(meta), path("*"), emit: dipspades_assembly
    path "versions.yml", emit : versions

    script:
    def log_file = "${meta.id}_dispades.log"

    def short_reads_param = "-1 ${R1} -2 ${R2}"
    
    script:
    """
    dipspades.py -v > versions.yml

    # Haplotype assembly - short reads -diploid genome - Illumina
    dipspades.py ${short_reads_param} \\
                ${extraOption} \\
                -o ${meta.id} > ${log_file} 2>&1
    """

}

