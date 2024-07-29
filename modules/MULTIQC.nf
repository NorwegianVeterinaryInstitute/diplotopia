//  Collecting results for report
process MULTIQC {
    conda "/cluster/projects/nn9305k/src/miniconda/envs/bifrost2022-multiqc"
    container 'quay.io/biocontainers/multiqc:1.20--pyhdfd78af_1'
    
    tag "$ID" 
    
    input:
    path(all)
    
    output: 
    file("*")
    
    script:
    """
    multiqc --version > multiqc.version
    multiqc --fullnames --config $projectDir/bin/multiqc_config.yaml $all 
    """

}
