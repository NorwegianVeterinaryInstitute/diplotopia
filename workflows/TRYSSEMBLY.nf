// workflow for multi-assembly testing
// make list assemblies 
// TODO - need to adjust inputs osv and params 


// modules 
//include { MASURCA } from "../modules/MASURCA.nf"
//include { NECAT_GLOBAL } from "../modules/NECAT.nf"
include {FOO_PATH} from "../modules/FOO.nf"


workflow TRYSSEMBLY {

    // here the pb is to standardize input 
    // so we need template and then the input is the location of those filled templates 

    input_ch = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header:['assembler', 'csv_input', 'comment'], skip: 1, sep:",", strip:true)
        .map { row -> (assembler, csv_input) =  [ row.assembler, row.csv_input ]}
    //[necat, /cluster/projects/nn9305k/active/evezeyl/projects/Saprolegnia/git/2023_Saprolegnia_pilot/dev_pipeline/tryssembly/input.necat.csv]
    //[masurca, /cluster/projects/nn9305k/active/evezeyl/projects/Saprolegnia/git/2023_Saprolegnia_pilot/dev_pipeline/tryssembly/input.masurca.csv    ]
    
    input_ch
    //.filter { row -> row.assembler == "/necat/" }
    .filter { assembler -> assembler == "/necat/" }
    .view()
    /*
    path_ch_necat = 
        input_ch.filter{ row.assembler == "necat" }
        .map { (necat_input_csv) = [it[1]] }

    FOO_PATH(necat_input_csv)
    FOO_PATH.out.necat_input_ch.view()
    */
        

        

    //masurca_input_ch = input_ch.filter{ row.assembler == "masurca" }

    
   
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


    
    // MASURCA 
    if (!params.assembly_list Contains masurca ) { 

        if (!params.input) { exit 1, "Missing input file"}

        // could modify to input directories containing short and long reads and assemblyID and concatenate
        input_ch = Channel
            .fromPath(params.input, checkIfExists: true)
            .splitCsv(header:['assemblyID', 'R1', 'R2', 'long_read'], skip: 1, sep:",", strip:true)
            .map { row -> (assemblyID, R1, R2, long_read) =  [ row.assemblyID, row.R1, row.R2, row.long_read ]}


        MASURCA_ASSEMBLY(input_ch)  
    }
    */ 



    

    

}
  

    


