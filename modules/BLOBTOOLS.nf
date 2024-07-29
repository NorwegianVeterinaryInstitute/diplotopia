// blobology https://github.com/DRL/blobtools
// https://blobtools.readme.io/docs/blobplot

process BLOBTOOLS {
    //conda "/cluster/projects/nn9305k/src/miniconda/envs/..."
    container 'quay.io/biocontainers/blobtools:1.1.1--py_1'
    
    //label 'process_high'
    
    input:
    path(assemblies)
    
    output: 
    path("quast_ref/*")


    script:
    """
    """
}



