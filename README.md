# NERVE

Group project with [MOLBINFO](http://www.bio.unipd.it/molbinfo/). 
[NERVE](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1570458/) means '*New Enhanced Reverse Vaccinology Environment*', and the purpose of this project is to update it and develop new modules. In particular this project includes:
- Perl to Python translation
- obsolete programs substitution
- [SPAAN](https://github.com/nicolagulmini/spaan) C to Python translation and improvement
- others... 

## REQUIREMENTS
NERVE accepts a prokaryotic proteome of a bacterium of which the user has to know if it is gram positive or gram negative.
Before starting to use NERVE, you should check if you have all the dependencies. Being a Python program, it imports some libraries, such as:
- Pandas 
- Bio (Biopython)
- Tensorflow

Also, you should do:
- `apt-get install ncbi-blast+` for blastp comparisons;
- `git clone https://github.com/Superzchen/iFeature` for features computation;
- `git clone git://github.com/nicolagulmini/NERVE` to import the modules (could be useful to periodically remove and reinstall this folder, to keep the program updated. To remove it, it is sufficient to run `rm -r NERVE` before cloning it);
- `python -m pip install git+https://github.com/nicolagulmini/tmhmm.py` for the third module, which needs tmhmm to compute the transmembrane domains. 

## Usage

Before describing each module, here a brief description of the parameters to pass. Some of them will be clear during the in detail description of the single modules:
- `-proteome1` (mandatory): the path to `.fasta` proteome file;
- `-gram` (mandatory): the gram (`p` for positive, `n` for negative or `a` for archea) of the organism.

Then the following parameters are optional:
- `-proteome2` (default = `None`) ...
- `-gram2` (default = `None`) ...
- `-psortb_output_path` (default = `None`) ...
- `-p_ad_no_citoplasm_filter` (default = `0.46`) ...
- `-p_ad_extracellular_filter` (default = `0.38`) ...
- `-transmemb_doms_limit` (default = `3`) ...
- `-percentage_of_covered_protein_for_razor` (default = `0.9`) ...
- `-e_value` (default = `1e-10`) ...
- `-similarity_function` (default = `0.8`) ...
- `-verbose` (default = `0`). Set to `1` (or any other symbol different from zero) if you what the program to print the protein information during the computation. Pay attention: could be a lot of data!

Note that some of the listed parameters are involved to the final scoring of the proteins, so changing them could deeply change the output of the program: be sure of what you do!

## Module 1: Subcelloc (still do not know the time to perform the computation on a realistic proteome)
For this module you will need the only external dependency of NERVE (that we are planning to substitute with a more convenient solution): [PSORTB](https://www.psort.org/psortb/). 
Here a very brief tutorial to have a working PSORTB from Linux command line, through Docker. The sequence of commands, starting without Docker, is:

```
sudo snap install docker
sudo docker pull brinkmanlab/psortb_commandline:1.0.2
wget https://raw.githubusercontent.com/brinkmanlab/psortb_commandline_docker/master/psortb
chmod +x psortb
```

once you have PSORTB, you are able to produce a file with the subcellular localization prediction of each of the proteins in your `.fasta` input file. If you do it before launching NERVE, then you should pass to it the path to the output file, specifying it to the `-psortb_output_path` parameter. Note: **you must produce a `terse` output**: if your input file is `input.fasta`, of a gram positive bacteria, in the psortb folder open the terminal and type
```
./psortb -i input.fasta -p -r ./ -o terse
```
where if your proteome is of a gram negative or archea bacteria, instead of putting `p` you should put `n` or `a`; and `-r ./` tells psortb to produce the output file in the same folder. 

## Module 2: Adhesin (approx 8 min for a realistic proteome)
You can find more details about this module [here](https://github.com/nicolagulmini/spaan). Here, for each protein, a probability value is computed, through an already trained neural network. That probability is the probability of the protein to be an adhesin. Take into account that the neural network has an accuracy of about 80%. Further improvements will be made, and this is one of the reasons why you should keep updated the nerve folder. 

## Module 3: Tmhelices (approx 4 min for a realistic proteome)
In this module are computed the transmembrane domains, thanks to the python library [tmhmm.py](https://github.com/dansondergaard/tmhmm.py).

### Module 3.1: Razor

## Module 4: Autoimmunity

### Module 4.1: Mouse immunity

## Module 5: Conservation (optional)

## Module 6: Function

## Module 7: Virulence
You can find more details about this module [here](https://github.com/nicolagulmini/virulent_factor_classification). In this case, for each protein some features are computed through `iFeature`, and then given to the already trained neural network in order to compute the probability to be a virulence factor. 

## Module 8: Select
This is the final module: if other modules will be added, they will be added before this module, because the final ranking of the proteins must take into account all the computed information. 
