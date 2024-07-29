// (population) Variants At Regions Where Reference Reads Uniquely Map module
include { MAPPING_VAR } from '../subworkflows/MAPPING_VAR.nf'
include { FREEBAYES } from '../subworkflows/FREEBAYES.nf'
include { TRIM } from '../modules/TRIM.nf'


workflow VARWRRUM {

    // input ch has a bit different format - now only illumina

    input_ch = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header:['ID', 'R1', 'R2', 'ploidy', 'comments'], skip: 1, sep:",", strip:true)
        .map { row -> (ID, R1, R2, ploidy) =  [row.ID, row.R1, row.R2 , row.ploidy]}

    illumina_input_ch = 
        input_ch
        .map { (ID, R1, R2) =  [it [0], it[1], it[2]] }

   
    ref_ch = Channel.fromPath(params.haploref, checkIfExists: true)




    // ---------------- OPTION using HaploSSembly -----------------------
    // reference given in parameter - assumed happlotype (! repeats non masked - maybe should)

    // Cleaning reads subworfklow
    TRIM(illumina_input_ch)

    // Illumina mapping subworflow 
    
    MAPPING_VAR(TRIM.out.reads, ref_ch) 

    // variant calling subworkflow
    //TODO if does not use ploidy should be removed 
    bam_ch = 
        input_ch
            .map { (ID, ploidy) = [it[0], it[3]]}
            .combine( MAPPING_VAR.out.bam, by: 0) 


    FREEBAYES(bam_ch, MAPPING_VAR.out.bam_bai, ref_ch) // channel [val(ID), val(ploidy), path(vcf)] AND ref_ch is path to fasta 
    /*
    
  
    - SEE calling variants in a population : how to do that 


    run freebayes - adapt options for diploid only
    https://bioinformaticschool.com/variant-calling-in-a-population-with-freebayes-maximizing-discriminant-power/


    output filtering: vcffilter -f "QUAL > 20" >results.vcf
    https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionVI/lessons/02_variant-calling.html
    
    QUAL and Depth DP or observation count (AO)
    QUAL : probability that there is a polymorphism at the loci described by the record
    1 - P(locus is homozygous given the data). 


    // ----------------OPTION using normal reference ---------------------
    // where we can get the one with coverage osv eventually on non called sites

    */



}