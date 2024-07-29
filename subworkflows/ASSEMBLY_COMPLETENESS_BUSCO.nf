// ASSEMBLY_COMPLETENESS_BUSCO.nf

include { BUSCO } from "../modules/BUSCO.nf"
include { BUSCO_REPORT } from "../modules/REPORTS.nf"

workflow ASSEMBLY_COMPLETENESS_BUSCO {
    take: 
    lineages // channel val id_lineage, path_download: list of different lineages and associated path
    assembly // channel val ID, path assembly 
    
    
    main:
    // Create combination : 
    // BUSCO will run for each lineage provided in params - for each assembly
    // tuple (id_assembly, path_assembly, id_lineage, path_download)
    inbusco = assembly.combine(lineages)
 
    BUSCO(inbusco)

    // getting only all the json files for summary
    all_json = BUSCO.out.busco_ch
                .map { paths = it[2] }
                .flatten()    
                .filter( ~/.*.json/ ) 
                .collect()

    BUSCO_REPORT(all_json)
        
    //emit:
}