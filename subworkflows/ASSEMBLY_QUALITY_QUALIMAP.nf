// Assembly quality evaluation workflow - with qualimap
include { MINIMAP } from "../modules/MINIMAP.nf"
include { SAMTOOLS_TOBAM_LONG; SAMTOOLS_MERGE_LONG_SHORT; SAMTOOLS_INDEX_SHORT } from "../modules/SAMTOOLS.nf"
include { BWAMEM } from "../modules/BWAMEM.nf"
include { QUALIMAP; QUALIMAP as QUALIMAPLONG; QUALIMAP as QUALIMAPSHORT } from "../modules/QUALIMAP.nf"
include { MULTIQC } from "../modules/MULTIQC.nf"

workflow ASSEMBLY_QUALITY_QUALIMAP {
    take:
    input_ch // chanel ID, assembly, longreads, R1, R2, type) 
    

    main: 
    // preparation: splitting for long reads mapping
    input_long_ch = input_ch
            .filter {it[5] == "long" | it[5] == "hybrid"}
            .map { (ID, assembly, longreads) =  [ it[0], it[1], it[2] ]}

    // splitting for short reads mapping
    input_short_ch = input_ch
            .filter {it[5] == "short" | it[5] == "hybrid"}
            .map { (ID, assembly, R1, R2) = [ it[0], it[1], it[3], it[4] ]}   



    // 1. MAPPING - Mapping of long and short reads, and hybrid assembly
    // 1.1 long reads
    MINIMAP(input_long_ch)
    SAMTOOLS_TOBAM_LONG(MINIMAP.out.longsam_ch)
            // note id for the different qualimaps need to be slighly different 
            // otherwise multiqc will not report all different names for hybrid, short or long

    // 1.2 short reads
    BWAMEM(input_short_ch)
    SAMTOOLS_INDEX_SHORT(BWAMEM.out.shortbam_ch)

    // 1.3 hybrid assembly - combining mapping short and long reads, using ID
    combinedbam_ch = 
        SAMTOOLS_TOBAM_LONG.out.longbam_ch
        .combine(SAMTOOLS_INDEX_SHORT.out.shortindexbam_ch, by: 0)
    
    SAMTOOLS_MERGE_LONG_SHORT(combinedbam_ch)

    // 2. QUALITY ASSEMBLY - Coverage statistics - per read type
    // I want part info provided both by short and long reads in case of hybrid assembly 
    // to have idea coverage individually (so qualimap done of each and combination)
        
    // 2.1 long reads
    QUALIMAPLONG(
        SAMTOOLS_TOBAM_LONG.out.longbam_ch.map { (ID, path) = [it[0] + "_long", it[1]] }
        ) 
 
    // 2.2 short reads
    QUALIMAPSHORT(
        SAMTOOLS_INDEX_SHORT.out.shortindexbam_ch
        .map { (ID, path) = [it[0] + "_short", it[1]] }
        )
    
    // 2.3 hybrid assembly
    QUALIMAP(
        SAMTOOLS_MERGE_LONG_SHORT.out.mergedbam_ch
        .map { (ID, path) = [it[0] + "_hybrid", it[1]] }
    
    ) 
    // 2.4 Report - Combining results for all qualimaps
    multiqc_ch = 
        QUALIMAPLONG.out.qualimap_ch.collect()
        .merge(QUALIMAPSHORT.out.qualimap_ch.collect())
        .merge(QUALIMAP.out.qualimap_ch.collect())


    MULTIQC(multiqc_ch)
    //emit:

}