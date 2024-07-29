include { QUAST_ASSEMBLY; QUAST_REF } from "../modules/QUAST.nf"

workflow ASSEMBLY_QUALITY_QUAST {
    take: 
    assemblies // chanel path assemblies - all emitted at once
    
    main:
    // Conditional quast type - with external reference: ref of own assembly as ref : assembly, or both
    if ( (params.which_quast=="both" | params.which_quast=="assembly") ) {
        
        // ensure necessary params are provided
        if ( params.quast_ref_assembly==null)  
        exit 1, "Please review your options:\
        Provide the path of the assembly that will be used as reference for QUAST."

        QUAST_ASSEMBLY(assemblies)
    }

    if ( (params.which_quast=="both" | params.which_quast=="ref") ) {
        
        // ensure necessary params are provided
        if ( params.quast_ref==null  | params.quast_annot==null)  
        exit 1, "Please review your options:\
        Provide an external reference AND its annotation,\
        or review if you are using the correct option to run QUAST."

        
        QUAST_REF(assemblies)

    }

//emit:
}