# NERVE

Group project with [MOLBINFO](http://www.bio.unipd.it/molbinfo/). 
[NERVE](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1570458/) means 'New Enhanced Reverse Vaccinology Environment'. The purpose of this project is to update it and to develop new modules. In particular this project includes:
- Perl to Python translation
- obsolete programs substitution
- [SPAAN](https://github.com/nicolagulmini/spaan) improvement
- others... 

# Pipeline

NERVE accepts a prokaryotic proteome of a bacterium of which the user has to know if it is gram positive or gram negative.

## Module 1: Subcelloc
For this module you will need the only external dependency of NERVE (that we are planning to substitute with a more convenient solution): PSORTB. 
Here a very brief tutorial to have a working PSORTB from Linux command line, through Docker. The sequence of commands, starting without Docker, is:

```
sudo snap install docker
sudo docker pull brinkmanlab/psortb_commandline:1.0.2
wget https://raw.githubusercontent.com/brinkmanlab/psortb_commandline_docker/master/psortb
chmod +x psortb
```
once you have PSORTB, you are able to produce a file with the subcellular localization prediction of each of the proteins in you `.fasta` input file. If you do it before launching NERVE, then you should pass to it the path to the output file. Note: you must produce a `terse` output!

## Module 2: Adhesin
## Module 3: Tmhelices
### Module 3.1: Razor
## Module 4: Autoimmunity
### Module 4.1: Mouse immunity
## Module 5: Conservation (optional)
## Module 6: Function
## Module 7: Virulence factor
## Module 8: Select
