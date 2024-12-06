# README 


using quarto in container

https://github.com/quarto-dev/quarto-cli/pkgs/container/quarto
```bash
cd $USERWORK/images
apptainer pull quarto_172.sif docker://ghcr.io/quarto-dev/quarto:1.7.2

IMG="/cluster/work/users/evezeyl/images/quarto_172.sif"
# this starts the container you will see change of prompt then
apptainer shell $IMG
```
