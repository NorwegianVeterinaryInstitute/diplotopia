// renaming assemblies to avoid name collision and allow to collect all assemblies for quast
process RENAME {

    input:
    tuple val(ID), path(assembly)
    
    output: 
    tuple val(ID), path("*.fasta"), emit: renamed_ch
    
    script:
    """
    mv $assembly ${ID}.fasta
    """
}