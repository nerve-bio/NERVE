# NERVE

Group project with [MOLBINFO](http://www.bio.unipd.it/molbinfo/). 
[NERVE](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1570458/) means '*New Enhanced Reverse Vaccinology Environment*', and the purpose of this project is to update it and develop new modules.

## REQUIREMENTS
NERVE accepts a prokaryotic proteome of a bacterium of which the user has to know if it is gram positive or gram negative.
Before starting to use NERVE, you should check if you have all the dependencies. Being a Python program, it imports some libraries, such as:
- Pandas 
- Numpy
- Bio (Biopython)
- Tensorflow (**pay attention: you must install the 2.6.0 version!**)

Also, you should do:
- `apt-get install ncbi-blast+` for blastp comparisons;
- `git clone https://github.com/Superzchen/iFeature`
- `git clone git://github.com/nicolagulmini/spaan` for features computation and trained models;
- `git clone git://github.com/nicolagulmini/NERVE` to import the modules (could be useful to periodically remove and reinstall this folder, to keep the program updated. To remove it, it is sufficient to run `rm -r NERVE` before cloning it);
- `python3 -m pip install git+https://github.com/nicolagulmini/tmhmm.py` for the third module, which needs tmhmm to compute the transmembrane domains. 

**Please be sure that all of these operations are done inside the same folder, in order to correctly import all the required libraries and modules.**

## Usage

Before describing each module, here a brief description of the parameters to pass. Some of them will be clear during the in detail description of the single modules:
- `-proteome1` (mandatory): the path to `.fasta` proteome file;
- `-gram` (mandatory): the gram (`p` for positive, `n` for negative or `a` for archea) of the organism.

Then the following parameters are optional:
- `-proteome2` (default = `None`) ...
- `-gram2` (default = `None`) ...
- `-razor` (default = `False`). Put `True` if you want the razor module to be called.
- `-razor_len` (default = `50`). If `razor` is `True`, then the `razor_len` is considered. In the razor module description the function of this parameter will be clarified.
- `-psortb_output_path` (default = `None`) ...
- `-p_ad_no_citoplasm_filter` (default = `0.46`) ...
- `-p_ad_extracellular_filter` (default = `0.38`) ...
- `-transmemb_doms_limit` (default = `3`) ...
- `-percentage_of_covered_protein_for_razor` (default = `0.9`) ...
- `-e_value` (default = `1e-10`) ...
- `-similarity_function` (default = `0.8`) ...
- `-verbose` (default = `0`). Set to `1` (or any other symbol different from zero) if you what the program to print the protein information during the computation. Pay attention: could be a lot of data!

Note that some of the listed parameters are involved to the final scoring of the proteins, so changing them could deeply change the output of the program: be sure of what you do!
Moreover it can also happen that the program asks you for the `[sudo] user password`, it does this to launch the commands that will be described later.

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

Also, we report the limitations of psortb:
- *Proteins resident at multiple localization sites*. Many proteins can exist at multiple localization sites. Examples of such proteins include integral membrane proteins with large periplasmic domains, or autotransporters, which contain an outer membrane pore domain and a cleaved extracellular domain. The current version of PSORTb handles this situation by flagging proteins which show a distribution of localization scores favouring two sites, rather than one. It is important to examine the distribution of localization scores carefully in order to determine if your submitted protein may have multiple localization sites and if so, which two sites are involved.
- *Lipoproteins*: The current version of PSORTb does not detect lipoprotein motifs.
- *Precision vs. Recall*: PSORTb is designed to emphasize precision (or specificity) over recall (or sensitivity). Programs which make predictions at all costs often provide incorrect or incomplete results, which can be propagated through annotated databases, datasets and reports in the literature. We believe that a confident prediction is more valuable than any prediction, and we have designed the program to this end. Note, however, that a user may choose to use their own reduced cutoff score in generating final predictions.

## Module 2: Adhesin (approx 8 min for a realistic proteome)
You can find more details about this module [here](https://github.com/nicolagulmini/spaan). Here, for each protein, a probability value is computed, through an already trained neural network. That probability is the probability of the protein to be an adhesin. Take into account that the neural network has an accuracy of about 80%. Further improvements will be made, and this is one of the reasons why you should keep updated the nerve folder. 

## Module 3: Tmhelices (approx 4 min for a realistic proteome)
In this module are computed the transmembrane domains, thanks to the python library [tmhmm.py](https://github.com/dansondergaard/tmhmm.py). 

### Module 3.1: Razor
In this module only the proteins with at least 3 transmembrane domains are considered.
For outermembrane proteins consider both the 'i' and 'o' loops, otherwise only the 'o' loop. 
Then take the longest out-membrane piece and replace the original sequence with it to perform the following analyses
(only if the longest piece is reasonably long...).

## Module 4: Autoimmunity (approx 10 min for blastp on sapiens, 10 sec for parsing, >30 min for comparison

### Module 4.1: Mouse immunity

## Module 5: Conservation (optional)

## Module 6: Function
Still implementing it. We want to include [DeepGO](https://github.com/bio-ontology-research-group/deepgo) in this module.

## Module 7: Virulence
You can find more details about this module [here](https://github.com/nicolagulmini/virulent_factor_classification). In this case, for each protein some features are computed through `iFeature`, and then given to the already trained neural network in order to compute the probability to be a virulence factor. 

## Module 8: Select
This is the final module: if other modules will be added, they will be added before this module, because the final ranking of the proteins must take into account all the computed information. 
