// Assembly contamination detection workflow
// To be improved 

include { KRAKEN as KRAKEN_ASSEMBLY } from "../modules/KRAKEN.nf"

workflow DETECT_ASSEMBLY_CONTAMINANTS {
    take:
    assembly // channel: [tuple val(ID), path(assembly)] 

    main: 
    // Detection contamination contigs tweak kraken 
    KRAKEN_ASSEMBLY(assembly)

    // Rscript to add the sankey plots correcting by size assembly for visualisation

    //emit:

}