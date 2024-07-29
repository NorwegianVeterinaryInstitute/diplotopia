# Developpment testing 
# TRYSSEMBLY 


## NECAT

cd /cluster/projects/nn9305k/active/evezeyl/projects/Saprolegnia/git/2023_Saprolegnia_pilot/dev_pipeline/tryssembly



MAIN="/cluster/projects/nn9305k/active/evezeyl/projects/Saprolegnia/git/2023_Saprolegnia_pilot/code/nf/main.nf" 
NF="/cluster/projects/nn9305k/bin/nextflow_23.04.4"
CONFIG="/cluster/projects/nn9305k/active/evezeyl/projects/Saprolegnia/git/2023_Saprolegnia_pilot/code/nf/conf/saga_sapro.config"

#run with apptainer - did not do conta 
module purge 
module load Java/17.0.4
$NF run $MAIN -c $CONFIG --out_dir . -work-dir $USERWORK/tryssembly --track tryssembly --input input_tryssembly.csv  -profile apptainer -resume | tee nf.runlog 

