// variant filtering
// Thanks to Snippy 
// expressions https://samtools.github.io/bcftools/bcftools.html#expressions
process BCFTOOLS_FILTER {
    //conda "/cluster/projects/nn9305k/src/miniconda/envs/"
    container 'quay.io/biocontainers/bcftools:1.9--ha228f0b_4'
    
    label 'process_high_memory_time2'
    tag "$ID"

    input:
    tuple val(ID), path(unfiltered)

    output: 
    tuple val(ID), path("*filtered.vcf"), emit: vcf
    path("*")

    script: 
    oldtxt = 'ID=QR,Number=1,Type=Integer'
    newtxt = 'ID=QR,Number=1,Type=Float'
    """
    # need to fix headers
    sed -i "s/$oldtxt/$newtxt/g" $unfiltered
 
    bcftools view $params.filter_args $unfiltered > ${ID}_filtered.vcf
    """
}

process BCFTOOLS_INDEX {
    //conda "/cluster/projects/nn9305k/src/miniconda/envs/"
    container 'quay.io/biocontainers/bcftools:1.9--ha228f0b_4'
    
    label 'process_high_memory_time2'
    tag "$ID"

    input:
    tuple val(ID), path(vcf)

    output: 
    tuple val(ID), path("*gz"), path("*.tbi"), emit: vcfgz_index
    path("*")

    script: 
    // using tbi for igv visualisation
    """
    bcftools convert -Oz -o ${vcf}.gz $vcf
    bcftools index --tbi -f ${vcf}.gz
    """
}

// creates a multi-sample file 
// https://www.biocomputix.com/post/how-to-combine-merge-vcf-bcf-files-using-bcftools-merge
process BCFTOOLS_MERGE  {
    //conda "/cluster/projects/nn9305k/src/miniconda/envs/"
    container 'quay.io/biocontainers/bcftools:1.9--ha228f0b_4'


    //label 'process_high_cpu_short'
    tag "$IDS"

    input:
    val(IDS)
    path(vcfgzs)
    path(tbis)
    path(ref)
    
    output: 
    tuple val(IDS), path("merged_indiv_calls.vcf.gz"), path("merged_indiv_calls.vcf.gz.tbi"), emit: vcfgz

    script:
    ID = IDS.join("_") 
    """
    bcftools merge --threads $task.cpus --gvcf $ref \\
    -Oz -o merged_indiv_calls.vcf.gz *.vcf.gz 
    bcftools index --tbi -f merged_indiv_calls.vcf.gz

    echo "merging done with IDS $ID" > bcftools_merge.log
    """
}
