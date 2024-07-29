# Saprolegnia - images build 

## Necat 

```shell
docker buildx build -f "./necat_dockerfile" -t "necat" . 2>&1 | tee build.log
docker run -it necat bash

# perl NECAT/Linux-amd64/bin/necat.pl 

sudo docker login
docker tag necat evezeyl/necat
docker push evezeyl/necat

singularity pull necat.sif docker://evezeyl/necat
singularity shell necat.sif --bind "directory_to_bind"

```


```
docker pull ttubb/nanopolish
# Testing if can work on a non root container
singularity pull nanopolish.sif docker://ttubb/nanopolish
singularity shell nanopolish.sif --bind "directory_to_bind"
sudo singularity shell --writable nanopolish.sif

```


```shell
cwltool --singularity https://github.com/common-workflow-language/common-workflow-language/raw/main/v1.0/v1.0/cat3-tool-mediumcut.cwl https://github.com/common-workflow-language/common-workflow-language/raw/main/v1.0/v1.0/cat-job.json
```


### [](https://github.com/common-workflow-language/cwltool?tab=readme-ov-file#running-a-tool-or-workflow-from-remote-or-local-locations)[  
](https://github.com/common-workflow-language/cwltool?tab=readme-ov-file#id10)

## Quast we need the databases


```shell
docker buildx build -f "./quast_dockerfile" -t "quast" . 2>&1 | tee build.log
docker run -it quast bash
conda env list # check versions 

sudo docker login
docker tag quast evezeyl/quast
docker push evezeyl/quast

singularity pull quast.sif docker://evezeyl/quast
singularity shell quast.sif  
# test need to be done
python 
import matplotlib 
import job
#singularity shell quast.sif --bind "directory_to_bind"
singularity shell --cleanenv quast.sif 

```



# Cleaning 
```
docker system prune -a
singularity cache clean
apptainer cache clean
```


# Notes important  - to put in manuals 
#container #important #notes
- path: container first, path rest at the end (avoid interaction environment - which can make that some libraries are not found eg. do not use the python that is installed in the container and then does not find the python libraries that actually were installed in the container !)
- building within container - best to go to directory with workdir and put the command there 
- 