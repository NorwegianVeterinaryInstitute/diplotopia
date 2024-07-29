// Concatenate results taxonomy with Blast and  
// Filter the contaminated contigs from each assembly


process FILTER_BLAST_CONTIGS { 

    //conda (params.enable_conda ? 'bioconda::bwa=0.7.8' : null)
	container 'evezeyl/compassr:latest'

    tag "${ID}"

    label 'process_medium_memory'


    input:
    tuple val(ID), path(contig_taxo), path(assembly) 
    // contig_taxo is paths for all contigs taxonomy
    // assembly is the renamed assembly

    output: 
    tuple val(ID), path("*_filtered_positive.fasta"), emit: decontassembly_ch 
    // in case all is the positive filter
    path("*_filtered_negative.fasta")
    path("*.{png,csv,html,command.sh,version,log,rds}")

    script: 
    javamem = "${task.memory.toGiga()-4}G"
    """
    cp $baseDir/bin/contigs_taxo_overview_filter.qmd .

    quarto render contigs_taxo_overview_filter.qmd \\
    -P res_dir:"." \\
    -P save_dir:"." \\
    -P ID:"${ID}" \\
    -P assemblyfile:"${assembly}" \\
    -P taxonomyDB:"$params.ranked_taxo_file" \\
    -P positive_filter:"$params.positive_filter" \\
    -P evalue_min:$params.evalue_min \\
    -P perc_identity_min:$params.perc_identity_min 2>&1 | tee ${ID}_compassr.log


    # need to rename the quarto render
    mv contigs_taxo_overview_filter.html ${ID}_contigs_taxo_overview_filter.html

    # copy the .command.sh to output if need adjustement
    cp .command.sh ${ID}.command.sh

    """

}
