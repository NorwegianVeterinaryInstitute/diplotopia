nextflow.enable.dsl=2
/*
module purge 
module load Java/17.0.4
NF="/cluster/projects/nn9305k/bin/nextflow_23.04.4"
$NF run test.nf
*/



// test8
ID = "TEST"
def contig = new File( 'TEST.part_contig2.fasta' ) 
contigID = contig.getBaseName(1) - ".part_" - ID
println contigID
// this works 

/*
//test 7

test_ch = Channel.of ( ["TEST_ID", ["path1", "path2", "path3"]])

// test_ch.view()
// TEST_ID
// [path1, path2, path3]
test_ch.transpose().view()

*/ 




// test 1
/*
params.lineage_dataset   =  "alveolata_odb10,stramenopiles_odb10"
params.dataset_path = "/cluster/projects/nn9305k/active/evezeyl/projects/Saprolegnia/git/2023_Saprolegnia_pilot/busco_downloads/busco_downloads/lineages"

lineages = params.lineage_dataset?.split(',') as List
path_lineages = lineages.collect { params.dataset_path + "/" +"$it" } 

println path_lineages


// Does not seems to be the way to go
// path_lineages2 = lineages.collect { (dirname, path_lineage) =  [ "$it", params.dataset_path + "/" +"$it" ]} 
//println path_lineages2

//path_lineages = lineages.flatten { params.dataset_path +"/" +"$it" }
//.join(' ') THIS IS TO CREATE STRINGS ! 



// checkif exist trough error when several paths apparently
lineages_ch = Channel.fromPath(path_lineages, type:'dir', checkIfExists: true)
    .flatten()
    .map { (id_lineage, path_lineage) = [it.name.split('_')[0], it] } 
    

lineages_ch.view()

//lineages_ch.getBaseName().view()

*/


// test 2
/*
numbers = Channel.of(1, 2, 3)
words = Channel.of('hello', 'ciao')

numbers
    .combine(words)
    .view()
*/


// test 3

/*
a = Channel.fromList( ['a1', 'b1', 'c1'])
b = Channel.fromList( ['a2', 'b2', 'c2'])
c = Channel.fromList( ['a3', 'b3', 'c3'])


// this does not keep the list - flattening
// c.concat( b, a ).view()

/*
a 
.merge(b)
.merge(c)
.view()

// [a1, a2, a3]
// [b1, b2, b3]
// [c1, c2, c3]


a1 = Channel.fromList( [['A', 'Aa1', 'Ab1'],['A', 'Aax', 'Abx']])
a2 = Channel.fromList( ['A', 'Aa2', 'Ab2'])
b1 = Channel.fromList( ['B', 'Ba1', 'Bb1'])




//a1.collect().view()

//a1.collect()
//.merge(a2.collect())
//.merge(b1.collect())
//.view()

a1
  .map { path = [it[1], it[2]]  }
  .collect()
  .view()
*/ 


/*
// test 4
      //  https://nextflow-io.github.io/patterns/process-when-empty/
        quast_ref_ch = 
        ? Channel.fromPath(params.quast_ref, checkIfExists: true)
        : Channel.empty()

        quast_annot_ch = 
        ? Channel.fromPath(params.quast_annot, checkIfExists: true)
        : Channel.empty()

*/ 


/*
// test 5
//params.quast_ref = null
//params.ref_assembly = null

params.quast_ref = "sth"
params.ref_assembly = "sth"


if ( (params.quast_ref==null & params.ref_assembly==null) | (params.quast_ref!=null & params.ref_assembly!=null)) 
        exit 1, "Please review your options:\
        Choose either an external or an internal reference for QUAST, but not both nor neither."
*/

// Test 6
/*

params.inputs = ""
  reads_ch = params.inputs
    ? Channel.fromPath(params.inputs, checkIfExists:true)

//test_ch.view()
is_empty = reads_ch.isEmpty()

if (is_empty) { println "Channel is empty" }
*/