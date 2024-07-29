// from Håkon
process TRIM {
	conda (params.enable_conda ? 'bioconda::trim-galore=0.6.10' : null)
	container 'quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0'

        label 'process_high'
        tag "$ID" 

        input:
        tuple val(ID), path(R1), path(R2)

        output:
        tuple val(ID), path {"*val_1.fq.gz"}, path {"*val_2.fq.gz"}, emit: reads

        script:
        """
        trim_galore -j $task.cpus --gzip -o . --paired --quality $params.phred_score -e $params.error_rate --length $params.minlength $R1 $R2 &> ${ID}_trimgalore.log
        """
}