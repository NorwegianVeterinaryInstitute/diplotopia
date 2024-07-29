process BUSCO {
    // does not exist conda yet
    conda "/cluster/projects/nn9305k/src/miniconda/envs/busco"
    container 'ezlabgva/busco:v5.7.1_cv1'

    label 'process_high_cpu_time'

    tag "$ID"
    
    input:
    tuple val(id_assembly), path(path_assembly), val(id_lineage), path(path_download)
    
    output: 
    tuple val(id_assembly), val(id_lineage), path("${id_assembly}_${id_lineage}/*"), emit: busco_ch
    path("*")
    
    script:
    // for online - not tested 
    if (params.lineage_dataset==null) {
        exit 1, "Solution for Offline option is not implemented yet - please choose the lineage dataset & provide the path"
        """
        busco -v > busco.version
        busco --cpu $task.cpus --in $path_assembly --auto-lineage --out ${id_assembly}_${id_lineage} --mode genome 
        """
    }
    // need to use $id_lineage because its the simplinked path 
    // for nf
    else {
        """
        busco -v > busco.version

        busco --cpu $task.cpus --in $path_assembly --out ${id_assembly}_${id_lineage} \\
        --mode genome \\
        --lineage_dataset $id_lineage --download_path $path_download --offline
        """
     }
}