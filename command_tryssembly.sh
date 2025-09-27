module purge 
module load Java/21.0.2  


# For masurca test
cd /cluster/projects/nn9305k/active/330403_001-A4_Tapeworm/analyses/assembly_masurca/20250926


DIPLO="/cluster/projects/nn9305k/active/evezeyl/projects/DIPLOTOPIA/diplotopia"
MAIN="${DIPLO}/main.nf" 
NF="/cluster/projects/nn9305k/bin/nextflow_25.04.7"
NFCONFIG="${DIPLO}/nextflow.config"
CONFIG="${DIPLO}/SAGA/saga_flatworm.config"
INPUT="${DIPLO}/tmp_dev/input_tryssembly.csv"



$NF run $MAIN -c $NFCONFIG -c $CONFIG --out_dir 582-1 −workDir $USERWORK/diplotopia/TRYSSEMBLY --input $INPUT -profile apptainer,tryssembly -resume

# 2>&1 | tee 2025-09-26_nf.runlog