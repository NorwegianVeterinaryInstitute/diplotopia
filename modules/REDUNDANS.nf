// https://github.com/Gabaldonlab/redundans/tree/master#
// Contains platanus assembler 


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
    # maybe need to stage file but then pb when missing files ... and paths
    # path of redundans.py in container is in root

    /root/src/redundans/redundans.py --version > versions.yml

    # important to put long_reads param at the end because create as string, in on new line produces error
    /root/src/redundans/redundans.py -v ${short_reads_param} ${long_reads_param} \\
                --mem ${task.mem} --tmp . \\
                -o ${meta.id} \\
                -t ${task.cpus} > ${log_file} 2>&1
    
    """

    // Specific mount
    extraOptions '-v /host/path/to/data:/container/path/to/data:ro' 
    // Or, for a named volume:
    // extraOptions '--mount source=my_volume,target=/container/path/to/data'


    // -v /host/path/to/data:/container/path/to/data:ro: This example bind-mounts a host directory /host/path/to/data to /container/path/to/data inside the container, with read-only permissions (ro).
    // --mount source=my_volume,target=/container/path/to/data: This example mounts a named Docker volume my_volume to /container/path/to/data inside the container.


}

//
// # if you wish to process local files, you need to mount the volume with -v
// ## make sure you are in redundans repo directory (containing test/ directory)
// docker run -v `pwd`/test:/test:rw -it cgenomics/redundans:latest /root/src/redundans/redundans.py -v -i test/*.fq.gz -f test/contigs.fa -o test/run1
