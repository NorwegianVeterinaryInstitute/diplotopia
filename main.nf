// SAPRO TOOLS 
nextflow.enable.dsl=2

//include { PURELY_RAW } from "./workflows/PURELY_RAW.nf"
include { TRYSSEMBLY } from "./workflows/TRYSSEMBLY.nf"
include { FILTER_CONTIGS } from "./workflows/FILTER_CONTIGS.nf"
include { COMPASS } from "./workflows/COMPASS.nf"
include { HAPLOPURGE } from "./workflows/HAPLOPURGE.nf"
include { VARWRRUM } from "./workflows/VARWRRUM.nf"

workflow{ 
	/*
	if (params.track == "purelyraw") {
		// sub_workflow - detectomg contamination from reads --- 
		// NOT ready 
		PURELY_RAW()
	}
	*/

	if (params.track == "tryssembly") {
		// sub_workflow - detectomg contamination from reads --- 
		// NOT ready 
		TRYSSEMBLY()
	}

	if (params.track == "filter_contigs") {
		// COMP-ASS for comparison assemblies 
		FILTER_CONTIGS() 
	}

    if (params.track == "compass") {
		// COMP-ASS for comparison assemblies 
		COMPASS() 
	}
	
	if (params.track == "haplopurge"){
		// Getting HaploSSembly from the best assembly 
		HAPLOPURGE()
	}

	if (params.track == "varwrrum"){
		// (population) Variants At Regions Where Reference Reads Uniquely Map 
		VARWRRUM()
	}


}

workflow.onComplete {
	log.info "".center(74, "=")
	log.info "Pipeline Complete!".center(74)
	log.info "Output directory: $params.out_dir".center(74)
	log.info "Duration: $workflow.duration".center(74)
	log.info "Nextflow version: $workflow.nextflow.version".center(74)
	log.info "".center(74, "=")
}

workflow.onError {
	println "Pipeline execution stopped with the following message: ${workflow.errorMessage}".center(74, "=")
}
