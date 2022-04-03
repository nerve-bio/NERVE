# Documentation

At the beginning of the developement, NERVE was a modular program, which means that it was built in several different modules. Now, the modules are still present but the computation is performed in a sequential manner. Each module computes some features on the same given proteome (probability of a protein to be an adhesin, or comparisons among the proteins of another given proteome, for instance), with the possibility to decide which modules to activate. Only at the end of the entire computation, the proteins of the proteome are returned, sorted according to the computed score, indicating the best candidates to produce a vaccine.

## Modules

### Intro
This module installs all the dependencies of the program, such as [iFeature](https://github.com/Superzchen/iFeature), [tmhmm.py](https://github.com/dansondergaard/tmhmm.py), [Biopython](https://biopython.org/), [PycURL](http://pycurl.io/), [DeepFRI](https://github.com/flatironinstitute/DeepFRI) and many others...

### Parameters
This module contains all the parameters (mandatory or not) which could be modified by the user, like the input, the verbose, or more technical ones like thresholds for models' predictions. 

### Subcelloc
Predicts the subcellular localization of the pathogen's proteins using [CELLO](http://cello.life.nctu.edu.tw/) predictor.

### Adhesin
Predicts the probability to be an adhesin protein with [ESPAAN](https://github.com/nicolagulmini/spaan).

### Tmhelices
Predicts the topology of your pathogen's proteins, in particular considering the number of transmembrane domains.

### Razor
Considers only the proteins with at least 3 transmembrane domains. 
For outermembrane proteins considers both the 'i' and 'o' loops, otherwise only the 'o' loop. 
Takes the longest out-membrane piece and replaces the original sequence with it to perform the following analyses
(only if the longest piece is reasonably long, let's say minimum 50 aa length)

### Autoimmunity 
Compares the pathogen's proteome with human proteome using BLASTp.Then, similar shared sequences are scanned in order to find MHC-I and MHC-II ligands. This is useful to help avoid possible autoimmunity reactions.

### Mouse immunity
Compares the pathogen's proteome with mouse proteome using BLASTp. Then, similar shared sequences are scanned in order to find MHC mouse ligands. This is useful to help avoid possible autoimmunity reactions in mice.

### Conservation
Compares the pathogen's proteome with another one from a different serogroup/strain, using BLASTp. The more the protein sequences are conserved among strains, the more protective the vaccine results.

### Virulence
Predicts the probability to be considered as a virulence factor, based on the pathogen protein sequences.

### Annotation
Functionally annotates the pathogen protein sequences retrieving info from Uniprot KB or using DeepFRI predictor, which assigns a GO term to each protein.

### Select
The final module which selects the proteine which are considered good candidates, and sorts them according to the similarity with the other given proteome. At the end, a .csv file is returned with all the computed information.
