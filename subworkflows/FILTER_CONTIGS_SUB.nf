// Sub-worflow to determine taxonomy of contigs in assemblies with blast
// and filter out contaminant contigs
// returns contaminants/host cleaned assembly
// https://sateeshperi.github.io/nextflow_varcal/nextflow/nextflow_sub_workflows


include { SEQKIT_TOCONTIG } from '../modules/SEQKIT.nf'
include { BLAST_CONTIG } from '../modules/BLAST_CONTIG.nf'
include { FILTER_BLAST_CONTIGS } from '../modules/FILTER_BLAST_CONTIGS.nf'

workflow FILTER_CONTIGS_SUB {

    take: 
    assembly // channel: [tuple val(ID), path(assembly)] "Channel where assemblies have been renamed with assembly ID"


    main: 
    // ------------- SPLIT assembly into CONTIGS -----
    SEQKIT_TOCONTIG(assembly)

    // ------------- BLASTS each contig to determine its taxonomy -----
    BLAST_CONTIG(
        SEQKIT_TOCONTIG.out.contigs_ch
        .transpose()
     ) 
    // BLAST_CONTIG.out.contigstaxo_ch // channel: [tuple val(ID), val(contigID), path("*taxo.tsv")] 
    

    // group taxo results for all contigs for each ID  
    blastout_ch = 
        BLAST_CONTIG.out.contigstaxo_ch
        .groupTuple(by: 0)
        .map { (ID, contigs) = [it[0], it[2]]}
        .combine(assembly, by : 0)
    // channel: [tuple val(ID), [path ID_contig_1, ..., path ID_contig_n], path(ID_assembly)] 
        
  
    // combine all results into one - filter out contaminants contigs - Rscript 
    FILTER_BLAST_CONTIGS(blastout_ch)
    
     
    
    // out : channel assembly filtered from contaminants
    emit: 
    // blast = BLAST_CONTIG.out.contigstaxo_ch // channel: [tuple val(ID), val(contigID), path("*taxo.tsv")]
    fasta = FILTER_BLAST_CONTIGS.out.decontassembly_ch // channel: [tuple val(ID), path(assembly_ok)] 
    

    
    }
    
