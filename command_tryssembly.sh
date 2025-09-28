# https://github.com/Gabaldonlab/redundans


tmux login4
 



# For masurca test
tmux
cd /cluster/projects/nn9305k/active/330403_001-A4_Tapeworm/analyses

module purge
module load Java/21.0.2 
DIPLO="/cluster/projects/nn9305k/active/evezeyl/projects/DIPLOTOPIA/diplotopia"
MAIN="${DIPLO}/main.nf" 
NF="/cluster/projects/nn9305k/bin/nextflow_25.04.7"
NFCONFIG="${DIPLO}/nextflow.config"
CONFIG="${DIPLO}/SAGA/saga_flatworm.config"
INPUT="${DIPLO}/tmp_dev/input_tryssembly.csv"
WORKDIR="/cluster/work/users/evezeyl/DIPLOTOPIA/TRYSSEMBLY"


$NF run $MAIN -c $NFCONFIG -c $CONFIG --out_dir diplotopia_assembly -work-dir $WORKDIR --input $INPUT -profile apptainer,tryssembly -resume

# 2>&1 | tee 2025-09-26_nf.runlog

srun --account=nn9305k --ntasks=1 --mem-per-cpu=4G --qos=devel --time=0:10:00 --pty bash -i
IMG=/cluster/work/users/evezeyl/images/quay.io-biocontainers-masurca-4.1.1--pl5321hb5bd705_0.img
apptainer shell $IMG
apptainer exec $IMG bash
find /usr/bin -name "masurca"
find /usr/local/bin -name "masurca" # /usr/local/bin/masurca
find /opt -name "masurca"

srun --account=nn9305k --ntasks=1 --mem-per-cpu=4G --qos=devel --time=0:10:00 --pty bash -i
module load MaSuRCA/4.1.0-GCC-11.3.0
# ok it seems the problem of running is when using the container 


IMG=/cluster/work/users/evezeyl/images/cgenomics-redundans-latest.img
apptainer shell $IMG
apptainer exec $IMG bash