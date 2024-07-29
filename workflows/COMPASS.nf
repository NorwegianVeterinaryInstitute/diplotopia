// Not in subworkflows
include { RENAME } from "../modules/RENAME.nf"

// Local: Sub-workflows
include { FILTER_CONTIGS_SUB } from '../subworkflows/FILTER_CONTIGS_SUB.nf'
include { DETECT_ASSEMBLY_CONTAMINANTS } from '../subworkflows/DETECT_ASSEMBLY_CONTAMINANTS.nf'
include { ASSEMBLY_QUALITY_QUALIMAP } from '../subworkflows/ASSEMBLY_QUALITY_QUALIMAP.nf'
include { ASSEMBLY_QUALITY_QUAST } from '../subworkflows/ASSEMBLY_QUALITY_QUAST.nf'
include { ASSEMBLY_COMPLETENESS_BUSCO } from '../subworkflows/ASSEMBLY_COMPLETENESS_BUSCO.nf'

workflow COMPASS {
    
    // ---------------------------------------------------------------------------------------------
    // ----------------------     INPUT:   preparation and manipulaton        ----------------------
    // ---------------------------------------------------------------------------------------------

    // input csv format ID, path assembly, path reads
    if (!params.input) { exit 1, "Missing input file"}

    input_ch = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header:['ID', 'assembly', 'longreads', 'R1', 'R2', 'type', 'comments'], skip: 1, sep:",", strip:true)
        .map { row -> (ID, assembly, longreads, R1, R2, type) =  [row.ID, row.assembly, row.longreads, row.R1, row.R2, row.type ]}
    
    // DEFENSIVE PROGRAMMING - INPUT CHECKS - TO ADD 

    // -------------------------- ASSEMBLIES ----------------------------------
    // renaming assemblies - avoiding name collision - emission all at same time

    assembly_ch = 
        input_ch 
        .map { (ID, assembly) =  [it[0], it[1]]}

    RENAME(assembly_ch)

    // ---------------------------------------------------------------------------------------------
    // ----------------------                    PIPELINE                     ----------------------
    // ---------------------------------------------------------------------------------------------
    

    
    
    // Defensive : option check 
    if ( params.remove_contaminants !="yes" && params.remove_contaminants !="no" ) {
        exit 1, """Please review your options:
        Determine if you only want to: 
        - option 'no':run directly the comparisson of assemblies hereby skipping contigs contamination removal
        (recommended if you did that before),
        OR
        - option 'yes': if you want to first remove contaminants and then compare assemblies
        (recommended at rerun only)."""
    }


    if ( params.remove_contaminants=="yes") { 

        FILTER_CONTIGS_SUB(RENAME.out.renamed_ch)

        assembly_ok_ch = FILTER_CONTIGS_SUB.out.fasta 

        } 
        else if (params.remove_contaminants=="no") {
            // out: bypassing - copy the raw assembly 
            assembly_ok_ch = RENAME.out.renamed_ch
        }

   
    // If only : then it should be finished - otherwise it continues to the rest of the pipeline
    if ( params.remove_contaminants=="yes" | params.remove_contaminants=="no" ) {
        
        DETECT_ASSEMBLY_CONTAMINANTS(assembly_ok_ch)

        // QUALIMAP - Assembly quality by mapping (decontaminated reads if choosen)  
           //replace the assembly that needs to be used if has been decontaminated 
           // channel: [(ID, assembly, longreads, R1, R2, type)] 
        new_input_ch = input_ch                 
                        .combine(assembly_ok_ch, by:0)
                        .map { (ID, assembly, longreads, R1, R2, type) =  [it[0], it[6], it[2], it[3], it[4], it[5]]}

        ASSEMBLY_QUALITY_QUALIMAP(new_input_ch)
        
      

        // QUALITY ASSEMBLY - Coverage statistics - QUAST
        quast_in_ch = assembly_ok_ch
                        .map { (assembly) = [ it[1] ] }
                        .collect()
        
        ASSEMBLY_QUALITY_QUAST(quast_in_ch)
    

        // QUALITY ASSEMBLY - COMPLETENESS - BUSCO 
        lineages_ch = 
            channel.from(params.lineage_dataset?.split(',') as List)
            .combine(Channel.fromPath(params.busco_download_path, type:'dir', checkIfExists: true))
            .map { (id_lineage, path_download) = [ it[0], it[1] ]}

        ASSEMBLY_COMPLETENESS_BUSCO(lineages_ch, assembly_ok_ch)
        

    // Could eventually be dotplot for ref and assemblies - for synteny - but need detail - dgenies is only interactive ?    
    }  


}
