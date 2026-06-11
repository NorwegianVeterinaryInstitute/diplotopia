#!/bin/bash
# Run filter_contigs track on SAGA
# Usage: bash SAGA/run_filter_contigs.sh

nextflow run main.nf \
    -c SAGA/saga_sapro.config \
    -profile apptainer \
    --track filter_contigs \
    --input samples.csv \
    --out_dir results \
    --blastDB /cluster/shared/databases/blast/nt \
    --ranked_taxo_file /path/to/ranked_taxo.tsv \
    --positive_filter "phylum == 'Oomycota'"
