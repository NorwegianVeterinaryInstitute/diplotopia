// WORFKLOW TO DETECT Contamination in reads

//include { BLOBTOOLS } from '../modules/BLOBTOOLS.nf'

//include { KRAKEN_READS  } from '../modules/KRAKEN_READS.nf'

//include {  } from '../modules/BLAST_CONTIG.nf'//
// include { READS_CONTAMINATION } from '../subworkflows/READS_CONTAMINATION.nf'
// see if 
    
    
    
workflow PURELY_RAW {

    if (!params.input) { exit 1, "Missing input file"}

    input_ch = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header:['ID', 'assembly', 'longreads', 'R1', 'R2', 'type', 'comments'], skip: 1, sep:",", strip:true)
        .map { row -> (ID, assembly, longreads, R1, R2, type) =  [row.ID, row.assembly, row.longreads, row.R1, row.R2, row.type ]}

    
    // Detection contamination from reads
    // KRAKEN_READS()  // eventually - add 

    
    // GC content vs Read coverage (short, long, hybrid if possible)
    // BLOBTOOLS()

    // remove contamination from reads 

    // ---- BLOBOLOGY on reads -----
    // GC content vs Read coverage (short, long, hybrid if possible)


}
