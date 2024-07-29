// variant calling with freebayes subworflow
// could try to paralellize leater on : look at  https://github.com/brwnj/freebayes-nf
// and improve the writing as nf core style ...

include { SAMBAMBA_MARKDUP } from "../modules/SAMBAMBA.nf"
include { SAMTOOLS_COV; SAMTOOLS_FAIDX; SAMTOOLS_INDEX_VAR } from "../modules/SAMTOOLS.nf"
include { FREEBAYES_CALL; FREEBAYES_REGION } from "../modules/FREEBAYES.nf"
include { BCFTOOLS_FILTER; BCFTOOLS_MERGE } from "../modules/BCFTOOLS.nf"
include { VT_NORMALIZE } from "../modules/VT.nf"

// for population 
include { FREEBAYES_CALL_POP } from "../modules/FREEBAYES.nf"
include { SAMTOOLS_MERGE_VAR; SAMTOOLS_SORT_VAR} from "../modules/SAMTOOLS.nf"

// variants index
include { BCFTOOLS_INDEX as INDEX_FREEBAYES_CALL;
          BCFTOOLS_INDEX as INDEX_BCFTOOLS_FILTER;
          BCFTOOLS_INDEX as INDEX_VT_NORMALIZE } from "../modules/BCFTOOLS.nf"

workflow FREEBAYES {
    take:
    bam // channel [val(ID), val(ploidy), path(bam)]
    fasta // channel: /path/to/ref.fasta  ref_ch

    
    main:
    //Part 1. mark duplicates (or remove duplicates)
    if (params.markdup == "mark" | params.markdup == "remove" ) {
        SAMBAMBA_MARKDUP(bam)
        //  tuple val(ID), val(ploidy), path("*_markdup.bam"), 
        new_bam = SAMBAMBA_MARKDUP.out.bam
        }
        else if (params.markdup == null) {
            new_bam = bam 
        } 
        else {
                exit 1, "Please review the parameters for marking or removing duplicates reads;\
                must be either 'mark' or 'remove' or 'null'."
        }
    // reindexing the bam files after marking duplicates per isolates
    SAMTOOLS_INDEX_VAR(new_bam) 
    // Part 2. Some statistics  (mapped on the reference)
    SAMTOOLS_COV(new_bam)  //SAMTOOLS_COV.out.coverage_stats // channel: [val(ID), path(coverage_data)]


    // Part 3. Variant calling
    // Indexing reference
    SAMTOOLS_FAIDX(fasta)
    // Getting region splits 
    FREEBAYES_REGION(SAMTOOLS_FAIDX.out.faidx_ch)


    // Variant calling per sample 
    FREEBAYES_CALL(
        SAMTOOLS_INDEX_VAR.out.bam
        .combine(SAMTOOLS_FAIDX.out.faidx_ch)
        .combine(FREEBAYES_REGION.out.regions)
        )

    // ---------------------------------------------------------------    
    // Variant calling population (only diploids for now)

    ids = SAMTOOLS_INDEX_VAR.out.bam 
        .map { (id) =  [it[0]] }
        .collect() 

    bams = SAMTOOLS_INDEX_VAR.out.bam 
        .map { (bam)  =  [it[2]] }
        .collect()

    bais = SAMTOOLS_INDEX_VAR.out.bam 
        .map { (bai)  =  [it[3]] }
        .collect()

    FREEBAYES_CALL_POP(ids, bams, bais, 
                       SAMTOOLS_FAIDX.out.faidx_ch,
                       FREEBAYES_REGION.out.regions)


    // ---------------------------------------------------------------    
    // merging bam file for all to be able to look at the population at once 
    // apparently doing that or above give same resuluts - merging to be able to look at one bam file in igv
    SAMTOOLS_MERGE_VAR(bams, bais)
    SAMTOOLS_SORT_VAR(SAMTOOLS_MERGE_VAR.out.bam) 


    // ---------------------------------------------------------------
    // Filter variants
    BCFTOOLS_FILTER(
        FREEBAYES_CALL.out.vcf
        .concat(FREEBAYES_CALL_POP.out.vcf)
    )

    // Normalize variants
    VT_NORMALIZE(BCFTOOLS_FILTER.out.vcf.combine(fasta))

    // ---------------------------------------------------------------
    // Indexing the different types of variant files (at the end) - for IGV visualisation 
    INDEX_FREEBAYES_CALL(
        FREEBAYES_CALL.out.vcf
        .concat(FREEBAYES_CALL_POP.out.vcf)
        )
    
    INDEX_BCFTOOLS_FILTER(BCFTOOLS_FILTER.out.vcf)
    INDEX_VT_NORMALIZE(VT_NORMALIZE.out.vcf) 


    // ---------------------------------------------------------------
    // Combine individual variant calls - for crosscheck with population call 
    // does not seems to work if uncompressed - so only do a the end 
    // here would have better have its by index rearranged but did not manage now
     
    vcfgzs = 
        SAMTOOLS_INDEX_VAR.out.bam
        .map { (id) =  [it[0]] }
        .combine(INDEX_VT_NORMALIZE.out.vcfgz, by: 0)

    vcfgzs_ids = vcfgzs 
        .map { (id) =  [it[0]] }
        .collect()

    vcfgzs_vcfgz = vcfgzs
        .map { (vcfgz)  =  [it[1]] }
        .collect()

    vcfgzs_index = vcfgzs
        .map { (index)  =  [it[2]] }
        .collect()

    BCFTOOLS_MERGE(vcfgzs_ids, vcfgzs_vcfgz, vcfgzs_index, fasta )                     
    
    // For now does not need to go further 
    //emit:
    //vcfgz = INDEX_VT_NORMALIZE.out. [val(ID), path(vcf.gz), path(*.tbi)]
    // to do after nf core way 
    //freebayes_version == FREEBAYES_CALL.version // path *.version 
    }
    
    
 
    



