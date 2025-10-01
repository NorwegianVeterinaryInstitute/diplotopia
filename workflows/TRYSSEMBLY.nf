// workflow for multi-assembly testing
// TODO - need to adjust inputs osv and params 


// modules 
include { MASURCA } from "../modules/MASURCA.nf"
include { DIPSPADES } from "../modules/DIPSPADES.nf"
include { REDUNDANS } from "../modules/REDUNDANS.nf"
include { PLATANUS_ASSEMBLE } from "../modules/PLATANUS.nf"
//include { NECAT_GLOBAL } from "../modules/NECAT.nf"
//include {FOO_PATH} from "../modules/FOO.nf"


workflow TRYSSEMBLY {

    // here the pb is to standardize input 
    // so we need template and then the input is the location of those filled templates 

    // SECTION inputs
    input_ch = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header:['id', 'assembler', 'path_reads','R1','R2','long_read', 'genome_size', 'advanced_options', 'ploidy', 'comment'], skip: 1, sep:",")      
        .map { row ->
    
            // for convenience usage at this stage need to be strings - because of need of configuration files for software
            def R1 = file("${row.path_reads}/${row.R1}", checkIfExists: true)
            def R2 = file("${row.path_reads}/${row.R2}", checkIfExists: true)
            def genome_size = row.genome_size ?: ''
            
            // Determine if long reads are provided
            def has_longRead = !!row.long_read
            def longRead = has_longRead ?
                            file("${row.path_reads}/${row.long_read}", checkIfExists: true) :
                            file("${workflow.projectDir}/assets/dummy.fa")

            
            def extraOption = row.advanced_options ?: ''

            // !! is a shorthand for a boolean check on a non-null string !!null -> F !!'astring' -> T
            def meta = [
                id: row.id, 
                assembler : row.assembler,
                genome_size : genome_size,
                ploidy : row.ploidy,
                has_longRead : has_longRead, 
                has_extraOption : !!row.advanced_options
            ]                    
            tuple(meta, R1, R2, longRead, extraOption)
            }
    
    //input_ch.view()
    // !SECTION


    // SECTION branching channel for different assemblers
    input_ch.branch {meta, R1, R2, longRead, extraOption ->
        go_to_masurca: meta.assembler == "masurca" 
        go_to_dispades: meta.assembler == "dipspades"
        go_to_redundans: meta.assembler == "redundans"
        go_to_platanus: meta.assembler == "platanus"
        // go_to_platanus_allee: meta.assembler == "platanus_allee"
        }
        .set { branched_ch }

    //branched_ch.go_to_masurca.view()
    // !SECTION



    // SECTION Assemblers  
    MASURCA(branched_ch.go_to_masurca)

    DIPSPADES(branched_ch.go_to_dispades
            .map{ meta, R1, R2, longRead, extraOption  -> 
                tuple(meta, R1, R2, extraOption)
                } 
            )
    REDUNDANS(branched_ch.go_to_redundans
            .map{ meta, R1, R2, longRead, extraOption -> 
                tuple(meta, R1, R2, longRead)
                } 
            )
    // platanus section
    PLATANUS_ASSEMBLE(branched_ch.go_to_platanus
            .map{ meta, R1, R2, longRead, extraOption -> 
                tuple(meta, R1, R2, longRead)
                } 
            )
    

    // !SECTION


}  


    


/*
    path_ch_necat = 
        input_ch.filter{ row.assembler == "necat" }
        .map { (necat_input_csv) = [it[1]] }

    FOO_PATH(necat_input_csv)
    FOO_PATH.out.necat_input_ch.view()
    */
        

        


    
   
   /*   

    if (params.assemblers == "necat") {
        
    } else {
        error "No assemblers specified"
    }
    


    if (!params.assembly_list Contains necat  ) { 

        input_ch = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header:['assemblyID', 'necatconf', 'readlist', 'comment'], skip: 1, sep:",", strip:true)
        .map { row -> (assemblyID, necatconf, readlist) =  [ row.assemblyID, row.necatconf, row.readlist ]}

        //input_ch.view()
        // here we must have something that build param file if possible to use memory and cpus

        // if (params.duplex) 
        // in modules is weird and needs to be fixed
        // all at once
        NECAT_GLOBAL(input_ch)

        //NECAT_CORRECT(input_ch)
        //NECAT_ASSEMBLE(NECAT_CORRECT.out.necat_correct_ch)
        //NECAT_BRIDGE(NECAT_ASSEMBLE.out.necat_assemble_ch)


    }


    

    */ 



    

    