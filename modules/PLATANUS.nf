// Platanus 1.2.4 
// Platanus constructs consensus sequence of homologous chromosomes at first and tries to split into each haplotype sequence.

// trim platanus_trim -i pe.fofn -t 16
// platanus_internal_trim -i mp.fofn -t 16
// SECTION  Warning 
// reads must have been trimmed beforeHAND or add option
// !SECTION


process PLATANUS_ASSEMBLE {
    //conda "/cluster/projects/nn9305k/src/miniconda/envs/"
    //conda (params.enable_conda ? 'bioconda::chewbbaca=3.3.1' : null)
    
    container 'cmonjeau/platanus:latest'
    
    label 'process_high_memory_time3'

    tag "$meta.id" 
    
    input:
    tuple val(meta), 
        path(R1, stageAs: 'R1.fastq'), 
        path(R2, stageAs: 'R2.fastq'), 
        path(longRead, stageAs: 'longRead.fa') 

    output: 
    tuple val(meta),
        path("${meta.id}_out/${meta.id}_contig.fa"), 
        path("${meta.id}_out/${meta.id}_contigBubble.fa"), 
        path("${meta.id}_out/${meta.id}_kmerFrq.tsv"), 
        path(R1), 
        path(R2), emit: platanus_assembled
    path "versions.yml", emit : versions
    // publishing results
    // path("${meta.id}_out"), emit: platanus_assembly_publish



    script:
    // SECTION define script outputs
    def log_file = "${meta.id}_platanus.log"
    def out_dir="${meta.id}_out"
    def memory = task.memory.toGiga().toInteger() 

    // def long_reads_param = '' 
    // if (meta.has_longRead) {
    //     long_reads_param = "--longreads ${longRead}"
    // }

    // SECTION MASURCA execution 
    """
    ## FIXME add option trimming step 
    ## platanus --version > versions.yml

    # This is the cleaned reads - for testing unziped before
    # gunzip -c R1.fq.gz > R1.fastq
    # gunzip -c R2.fq.gz > R2.fastq

    Platanus assemble -o ${meta.id} -f R[12].fastq \\
        -k 28 -s 4 -n 0 -c 2 -tmp . \\
        -t ${task.cpus} -m ${memory} > ${log_file} 2>&1

    ## Grouping all results
    mkdir -p ${out_dir}
    mv ${meta.id}* ${out_dir}/
    """
}

// process PLATANUS_SCAFFOLD {
//     //conda "/cluster/projects/nn9305k/src/miniconda/envs/"
//     //conda (params.enable_conda ? 'bioconda::chewbbaca=3.3.1' : null)
    
//     container 'cmonjeau/platanus:latest'
    
//     label 'process_high_memory_time'

//     tag "$meta.id" 
    
//     input:
//     tuple val(meta), 
//         path(contig, stageAs: "platanus_contig.fa"), 
//         path(contigBubble, stageAs: "platanus_contigBubble.fa"), 
//         path(R1, stageAs: 'R1.fastq'), 
//         path(R2, stageAs: 'R2.fastq') 

//     output: 
//     tuple val(meta), 
//     path("${meta.id}_out/${meta.id}_scaffold.fa"), emit: platanus_scaffolded
//     path "versions.yml", emit : versions

//     script:
//     // SECTION define script outputs
//     def log_file = "${meta.id}_platanus_scaffold.log"
//     def out_dir="${meta.id}_out"
//     def memory = task.memory.toGiga().toInteger()

//     """
//     ## heu there is no options of memory here ...
//     Platanus scaffold \\
//         -o ${meta.id} \\
//         -c ${meta.id}_contig.fa \\
//         -b ${meta.id}_contigBubble.fa \\
//         -ip1 R1.fastq \\
//         -op1 R2.fastq \\
//         -n1 100 \\
//         -a1 150 \\
//         -d1 25 \\
//         -v 15 \\
//         -t ${task.cpus} >> ${log_file} 2>&1

//     """ 



//     Platanus gap_close \\
//      -o ${meta.id} \\
//      -c ${meta.id}_scaffold.fa \\
//      -f single_end_file.fastq \\ cannot implement this for now
//      -ip1 R1.fastq \\
//      -op1 R2.fastq \\
//      -s1 32 \\
//      -k 32 \\
//      -vo 32 \\ 
//      -vd 32 \\
//      -tmp . \\
//      -t ${task.cpus} >> ${log_file} 2>&1


// PLACE Holder - we can try later 
// https://cell-innovation.nig.ac.jp/maser/ToolsList/P000001849.html
// Platanus-allee tries to construct each haplotype sequence from the beginning and pair them as homologous chromosomes, 
// process PLATANUS_ALLELE {
// platanus_allee
// }