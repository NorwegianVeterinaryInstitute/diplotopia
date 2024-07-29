// masking low complexity regions in the references
// https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmask-guide/
// https://github.com/BioInfoTools/BBMap/blob/master/sh/bbmask.sh

process BBMASK  {
    //conda "/cluster/projects/nn9305k/src/miniconda/envs/"
    container 'quay.io/staphb/bbtools:39.06'

    tag "$ID" 
    label 'process_high'

    input: 
    tuple val(ID), path(assembly)
    path(sam) // must be optional file - can be a dummy
    
    output: 
    path("*_ref_masked.fasta"), emit: ref 
    tuple val(ID), path("*_ref_masked.fasta"), emit: ref_ID

    script:
    if (params.bbmask_maxcov != null) {
        """
        bbmask.sh threads=$task.cpus in=$assembly out=${ID}_ref_masked.fasta \\
        maskrepeats=t kr=5 minlen=40 mincount=4 \\
        masklowentropy=t ke=5 window=80 entropy=$params.bbmask_entropy \\
        mincov=10 maxcov=$params.bbmask_maxcov \\
        sam=$sam

        bbmask.sh --version > bbtools.version
        """
            
    } else {
        // masks low complexity regions only
       """
        bbmask.sh threads=$task.cpus in=$assembly out=${ID}_ref_masked.fasta \\
        maskrepeats=t kr=5 minlen=40 mincount=4 \\
        masklowentropy=t ke=5 window=80 entropy=$params.bbmask_entropy 

        bbmask.sh --version > bbtools.version
        """
       
    }
    
}