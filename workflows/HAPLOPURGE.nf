include { HAPLO_MINIMAP } from '../modules/MINIMAP.nf'
include { HAPLO_TOBAM } from '../modules/SAMTOOLS.nf'
include { HAPLOPURGE_PURGE_HIST; PURGE_COVERAGE; PURGE_HAPLOTIGS } from '../modules/PURGE_HAPLOTIGS.nf'

// Purging assemblies -> to haploid to be sure 
// Then we can rerun BUSCO and QUAST - so the first pipeline compare 
// https://bitbucket.org/mroachawri/purge_haplotigs/src/master/

// NOW ONLY FOR LONG READS 

workflow HAPLOPURGE {

    // ensure necessary params of purge phase are provided
    if ( params.purge_phase==null  | !( params.purge_phase=="hist" | params.purge_phase=="cov" | params.purge_phase!=="purge" | params.purge_phase=="all")) {
        exit 1, "Please review your options:\
        Provide a parametter for haplopurge phase."
    }
        
    input_ch = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header:['ID', 'assembly', 'longreads', 'R1', 'R2', 'type', 'haplo_low', 'haplo_mid', 'haplo_high', 'comments'], 
            skip: 1, sep:",", strip:true
            )
        .map { row -> (ID, assembly, longreads, R1, R2, type, haplo_low, haplo_mid, haplo_high) = 
            [row.ID, row.assembly, row.longreads, row.R1, row.R2, row.type, row.haplo_low, row.haplo_mid, row.haplo_high]}

    // DEFENSIVE PROGRAMMING - INPUT CHECKS - TO ADD (maybe here can find the way for the NA handling)   


    
    // -------------------------- CHANNELS INPUT ----------------------------------
    // from long reads
    long_input_ch = 
        input_ch
            .filter { it[5] == "long" | it[5] == "hybrid" }
            .map { (ID, assembly, longreads) =  [ it[0], it[1], it[2] ] }

    // for short reads

    // parameters for purge haplotigs
    haplo_param_ch = 
        input_ch
           .map { (ID, haplo_low, haplo_mid, haplo_high) =  [ it[0], it[6], it[7], it[8] ] }
           //.filter (it[6].isInteger() ) ? find a way if NA do not add 

    
    // -------------------------- PIPELINE ----------------------------------

    // for long reads

    HAPLO_MINIMAP(long_input_ch) // out channel: [tuple val(ID), path(assembly), path(longsam)]
    HAPLO_TOBAM(HAPLO_MINIMAP.out.longsam_ch) // out channel: [tuple val(ID), path(assembly), path(longbam)]

    // TODO for short reads TODO

    HAPLOPURGE_PURGE_HIST(HAPLO_TOBAM.out.longbam_ch)  // out channel: [tuple val(ID), path(assembly), path(gencov)]

    // here would be good to find a way to control if haplo params are not given - did not find way yet 

    if (params.purge_phase=="cov" |  params.purge_phase=="purge" | params.purge_phase=="all") { 
        // does from step 2  
        // channel [tuple val (ID), path(assembly), path(gencov), val(haplo_low), val(haplo_mid), val(haplo_high)]
        PURGE_COVERAGE(HAPLOPURGE_PURGE_HIST.out.haplohist_ch.combine(haplo_param_ch, by: 0))
    }
    
    if ( params.purge_phase=="purge" | params.purge_phase=="all") { 
        // does from step 3 or all 
        //channel [tuple val (ID), path(assembly), path(coverage_stats), path(bam), path(bamindex)]
        PURGE_HAPLOTIGS(PURGE_COVERAGE.out.haplopurged_ch.combine(HAPLO_TOBAM.out.mapindex_ch, by: 0) )
  

    // TODO - Add a step of validation of what has been filtered out (eg. junk) - identify what is its
    // BLAST (same taxonomy blast)

    // TODO - Add a busco completness step check on the purded assemblies (or separate ? is needed when want to look at it)

    // TODO - Add a summary of converage contig and filtered contigs ? find better way to croscheck info ()and completeness ?)


    // output file to launch the compass pipeline eventually - using the haplotype assemblies
    // add the paths osv to the launch

    // https://stackoverflow.com/questions/76590728/writing-values-from-a-channel-to-a-file
    // Channel [tuple val(ID), path(haplosembly, path(longreads), path(R1), path(R2), val(type)] - then need to create string

    // we want to retake the results from the result directory not from the NF_WORKDIR so need some modification 
    // https://github.com/nextflow-io/nextflow/discussions/3173 getting the paths 

 
    /*
    // function to get the absolute path of a file and return into a string

    def getItAbsolutePath (string_path) {
        File f = new File(string_path)
        full_string_path = f.getAbsolutePath()
        return full_string_path
        
    }
   
    
    //def rootpath = new File(".")
    //println rootpath.absolutePath
    //.toString()

    //println params.out_dir.getAbsolutePath()

    //.map {it -> tuple( it[0], (new File(it[1])).absolutePath)}
    */
    
    // first we need to get the path of the output not of the one in NF workdir 
    // problem here - only get the relative path - tried to fix above but how ?
    // each line needs to be a string with 

    /*

    my_values_ch =  
        PURGE_HAPLOTIGS.out.haplopurged_ch
        .map { it -> tuple( it[0], "${params.out_dir}/results/" + "happlopurge/02_HAPLOPURGE/HAPLOTIGS/" + it[1].getName() )}
        .combine(input_ch, by: 0)
        .map { (ID, assembly, longreads, R1, R2, type) = [it[0], it[1], it[3], it[4], it[5], it[6] ] }
        .map( it ->  it.join(","))
        .view()
    
    my_header_ch = Channel
        .of( "ID,assembly,longreads,R1,R2,type,comments" )

    my_header_ch
        .concat( my_values_ch )
        .collectFile( name: 'input_for_compass.csv', newLine: true, sort: false , storeDir: "${params.out_dir}/results")
    }
    */ 

 }
}
    
         
    