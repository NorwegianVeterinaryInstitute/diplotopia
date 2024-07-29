// preparation of reference for variant calling
// low complexity region masking, high coverage region masking ...

include { TRIM } from '../modules/TRIM.nf'
include { BWAMEM_VAR } from '../modules/BWAMEM.nf'

include { 
        SAMTOOLS_TOBAM; 
        SAMTOOLS_NAMESORT;
        SAMTOOLS_FIXMATE; 
        SAMTOOLS_COORDSORT_INDEX; 
        SAMTOOLS_COORDSORT_INDEX as SAMTOOLS_COORDSORT_INDEX2;
        SAMTOOLS_MARKDUP; 
        SAMTOOLS_FLAGSTAT; 
        SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT2;
        SAMTOOLS_BAMSTATS;
        SAMTOOLS_BAMSTATS as SAMTOOLS_BAMSTATS1;
        SAMTOOLS_BAMSTATS as SAMTOOLS_BAMSTATS2;
        SAMTOOLS_FAIDX; 
        SAMTOOLS_FILTER;
        SAMTOOLS_COVERAGE;
        SAMTOOLS_COVERAGE as SAMTOOLS_COVERAGE2;
          } from '../modules/SAMTOOLS.nf'

include { PURGE_HIST } from '../modules/PURGE_HAPLOTIGS.nf'

workflow VAR_PREP_BASE {

    take:
    assembly_reads_input_ch // channel: tuple val(ID) path(assembly), path(R1), path(R2)

    main:  
    // index for reference 
    SAMTOOLS_FAIDX(assembly_reads_input_ch.map { (ID, assembly) = [it[0], it[1]] })

    // prepare reads for mapping
    TRIM(assembly_reads_input_ch.map { (ID, R1, R2) = [it[0], it[2], it[3]] })

    // mapping
    BWAMEM_VAR(
        TRIM.out.reads
        .combine(assembly_reads_input_ch.map { (ID, assembly) = [it[0], it[1]]}, by : 0)
       )

    SAMTOOLS_TOBAM(BWAMEM_VAR.out.sam)
    
    // Coverage info before filtering
    SAMTOOLS_COORDSORT_INDEX(SAMTOOLS_TOBAM.out.bam)
    SAMTOOLS_COVERAGE(SAMTOOLS_COORDSORT_INDEX.out.bam_bai)
    SAMTOOLS_BAMSTATS(SAMTOOLS_COORDSORT_INDEX.out.bam_bai)


    // preparation tags for filtering 
    SAMTOOLS_NAMESORT(SAMTOOLS_TOBAM.out.bam)
    SAMTOOLS_FIXMATE(SAMTOOLS_NAMESORT.out.bam)
    SAMTOOLS_COORDSORT_INDEX2(SAMTOOLS_FIXMATE.out.bam)
    SAMTOOLS_MARKDUP(SAMTOOLS_COORDSORT_INDEX2.out.bam_bai)
    SAMTOOLS_FLAGSTAT(SAMTOOLS_MARKDUP.out.bam_bai)
    SAMTOOLS_BAMSTATS1(SAMTOOLS_MARKDUP.out.bam_bai) // replace flagstat ?  // rename to 2 would be better
    
    // Filtering reads with low mapping quality, and secondary alignments
    // We need to select one reference+index 
    // It must be the same path selected for resume (its always the same ref anyway)

    noid_ref_index = SAMTOOLS_FAIDX.out.noid_assembly_fai.take( 1 )

    SAMTOOLS_FILTER(
        SAMTOOLS_MARKDUP.out.bam_bai
        .combine(noid_ref_index)
    )

    
  
    // TOCHECH export reads and remap - seems flags were gone ? reflag
    // Rechecking filters
    SAMTOOLS_FLAGSTAT2(SAMTOOLS_FILTER.out.bam_bai)
    SAMTOOLS_BAMSTATS2(SAMTOOLS_FILTER.out.bam_bai) // replace flagstat ? 

    // hist coverage to be able to easy visualize haplo/diploid
    PURGE_HIST(SAMTOOLS_FILTER.out.bam_tohist)
    // and coverage info after filtering 
    SAMTOOLS_COVERAGE2(SAMTOOLS_FILTER.out.bam_bai)

    // remove duplicates - deprecated 
    //SAMBAMBA_MARKDUP(bam)

    emit:
    bam_bai = SAMTOOLS_FILTER.out.bam_bai // tuple: [ID, path(bam), path(bai)]
    ref_index = noid_ref_index // path(ref) path(index)
 

    
}