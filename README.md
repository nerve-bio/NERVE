# NERVE

Group project with [MOLBINFO](http://www.bio.unipd.it/molbinfo/). 
[NERVE](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1570458/) means 'New Enhanced Reverse Vaccinology Environment'. The purpose of this project is to update it and to develop new modules. In particular this project includes:
- Perl to Python translation
- obsolete programs substitution
- [SPAAN](https://github.com/nicolagulmini/spaan) improvement
- others... 

## REQUIREMENTS
NERVE accepts a prokaryotic proteome of a bacterium of which the user has to know if it is gram positive or gram negative.
Before starting to use NERVE, you should check if you have all the dependencies. Being a Python program, it imports some libraries, such as:
- Pandas 
- Bio (Biopython)
- Tensorflow

Also, you should do:
- `apt-get install ncbi-blast+` for blastp comparisons;
- `git clone git://github.com/nicolagulmini/NERVE` to import the modules (could be useful to periodically remove and reinstall this folder, to keep the program updated. To remove it, it is sufficient to run `rm -r NERVE` before cloning it);
- `python -m pip install git+https://github.com/nicolagulmini/tmhmm.py` for the third module, which needs tmhmm to compute the transmembrane domains. 

## Usage

Before describing each module, here a brief description of the parameters to pass:
- `-proteome1` (mandatory): the path to `.fasta` proteome file;
- `-gram` (mandatory): the gram (`p` for positive or `n` for negative) of the organism;
Then the following parameters are optional:
- `-proteome2` (default = `None`)
- `-gram2`(default = `None`)
- `-p_ad_no_citoplasm_filter` (default = `0.46`)
- `-p_ad_extracellular_filter` (default = `0.38`)
- `-transmemb_doms_limit` (default = `3`)
- `-percentage_of_covered_protein_for_razor` (default = `0.9`)
- `-e_value` (default = `1e-10`)
- `-similarity_function` (default = `0.8`)
- `-verbose` (default = `0`). Set to `1` if you what the program to print the protein information during the computation (pay attention: could be a lot of data!).

Note that some of the listed parameters are involved to the final scoring of the proteins, so changing them could deeply change the output of the program: be sure of what you do!

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
