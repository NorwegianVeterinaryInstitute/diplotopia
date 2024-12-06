# README 


# Quarto on SAGA 
## using quarto in container
> usage was limited when used with vscode. But the rendering function well 

[link image ](https://github.com/quarto-dev/quarto-cli/pkgs/container/quarto)

```bash
cd $USERWORK/images
apptainer pull quarto_172.sif docker://ghcr.io/quarto-dev/quarto:1.7.2

IMG="/cluster/work/users/evezeyl/images/quarto_172.sif"
# this starts the container you will see change of prompt then
apptainer shell $IMG
```
## installation quarto from tarbal in home directory (should put that in conda I think)

```bash
cd /cluster/home/evezeyl/software
# instructions here https://quarto.org/docs/download/tarball.html?version=1.6.39&idPrefix=download
mkdir quarto
tar -C  quarto -xvzf quarto-1.6.39-linux-amd64.tar.gz
ln -s ~/software/quarto/quarto-1.6.39/bin/quarto  ~/bin/quarto
quarto check
# editing bashrc and adding alias
alias quarto="~/bin/quarto"
# then make sure the path is configured on vscode for saga and workspace in the quarto extension (but not for user as different on pc)
```
Ok still things laking for the preview render (I do not know why, maybe because latex not detected osv) 