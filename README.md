 # NERVE

<div align="center"
  
  ![nerve](https://github.com/Andrea0097/NERVE/assets/113541183/cebc9a7c-daf5-4ffb-8575-dc5c06df13af)

</div>
  
<div align="left" 
  
  NERVE (New Enhanced Reverse Vaccinology Environment) is an open-source, reverse vaccinology environment, with which you can analyze bacterial proteomes in FASTA format to get the best protein vaccine candidates (PVCs) and their epitopes. 
  <p>
  The project was initially developed in Perl in 2006. You can find all the related information on the article at: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1570458
  </p>
  <p>
  Now, we are carrying on this project independently from the original developer and this is the Github page of the NERVE 2.0-stand-alone version. 
  </p>
  <p>
  Users can also leverage a user-friendly web-based 
  graphical interface (GUI) hosted on our dedicated server: https://nerve-bio.org.  This GUI empowers even novice users to operate NERVE with ease, requiring only the proteome FASTA file as input. It further enhances the experience by allowing visualization of analysis results using our 
  custom-built visualizer. Additionally, this tool facilitates the identification of epitopes (and their constituent parts) predicted to bind to multiple HLA alleles.
  </p>
</div>

## Key features:

 * **Machine Learning and Alignment Analysis**, using  ML methods for antigen analysis
 * **Flexible Usage**, customizing your settings
 * **Stand-alone Docker Version**, reducing dependencies for its installation
 * **Intuitive Web-based GUI**,
 * **Complete Python Rewrite**, enhancing its performance and maintainability.
 * **Ongoing Development**, continuously evolving with new features and improvements.
 * **Comprehensive Output**, including:
    * A CSV file summarizing all candidate proteins with scores and predicted characteristics (P_adhesin, P_antigen, location, etc.).
    * A separate CSV file listing discarded proteins.
    * Dedicated folders containing additional files (CSVs and PNGs) for each selected protein with predicted linear epitopes. 

**Benefits for Vaccine Research:**

* **Accessibility:** Researchers worldwide can access and use NERVE without barriers, accelerating vaccine development efforts globally. 
* **Transparency:** Allowing researchers to scrutinize the code, identify potential issues, and contribute to its refinement. 
* **Collaboration:** By enabling researchers to share, modify, and extend the pipeline, leading to innovative advancements. 
* **Community-Driven Development:** The open-source model encourages community participation, ensuring that NERVE remains aligned with the evolving needs of the vaccine research community. 

We welcome contributions from researchers interested in shaping the future of NERVE! Whether it's reporting bugs, suggesting improvements, or developing new features, your input is valuable.


## Installation procedure of stand-alone version:
NERVE can be used as a stand alone version taking advantage of [Docker](https://www.docker.com/) in Linux systems.

1) Install Docker following [these instructions](https://docs.docker.com/engine/install/) and [the post-installation procedure](https://docs.docker.com/engine/install/linux-postinstall/)
2) Clone the repository:
```
git clone https://github.com/nerve-bio/NERVE.git
```
4) Navigate to the correct folder:
```
cd NERVE
```
4) Run NERVE (the first time it will take a few minutes)
```
./NERVE.sh --help

```

This is the expected output:
```
usage: NERVE.py [-h] [-a] [-ev] -g [-ml] [-mm] [-m] [-mpsl] -p1 [-p2] [-pl] [-rz] [-rl] [-s] [-ss] [-tdl] [-vl] [-vir]
                [-ep] [-m1l] [-m2l] [-m1ovr] [-m2ovr] [-prt] [-wd] [-nd] [-id] [-dfd] [-epp]

NERVE arguments:
  -h, --help        show this help message and exit
  -a, --annotation  Activation (True) or deactivation (False) of annotation module. Uses DeepFri to retrieve protein functional onthologies (default: True)
  -ev, --e_value    Expect-value used in blastp for immunity modules (default: 1e-10)
  -g, --gram        Negative (n) or positive (p) gram stain of the pathogen of interest (default: None)
  -ml, --minlength  Minimal length required for shared peptides to be extracted in comparison analyses versus human and/or mouse (default: 9)
  -mm, --mismatch   Maximal number of not compatible substitutions allowed in shared peptides alignment windows of 'minlength' size in immunity modules (default: 1)
  -m, --mouse       Activation (True) or deactivation (False) of the mouse immunity module. This module compares proteome1 with mouse proteome and a further analysis of the
                        eventual shared peptides is carried out as in the autoimmunity module (default: False)
  -mpsl, --mouse_peptides_sum_limit 
                        Parameter calculated in mouse module and used by select module. Protein with 'sum of shared peptides of the i-protein with mouse proteins/number of aminoacids
                        of the i-protein' <= mouse_peptides_sum_limit and with absence of match mhc-I and Mhc-II mouse ligands are selected (default: 0.15)
  -p1, --proteome1  Path to proteome or Uniprot proteome ID (see: https://www.uniprot.org/proteomes/?query=&sort=score) (default: None)
  -p2, --proteome2  Path to proteome or Uniprot proteome ID (see: https://www.uniprot.org/proteomes/?query=&sort=score) (default: None)
  -pl, --padlimit   Set the probability of adhesin (pad) value cut-off for all proteins in select module. Thus, these proteins with a pad value < cut-
                        off are discarded (0.-1) (default: 0.5)
  -rz, --razor      Activation (True) or deactivation (False) of the loop-razor module. This module allows the recovery of protein vaccine candidates, with more than 2
                        transmembrane domains, that would otherwise be discarded in the select module. The longest loop with minimum len == 'razlen' aa will replace the original
                        protein sequence for following NERVE steps (default: False)
  -rl, --razlen     Set minimal length of loop considered in loop-razor module (default: 50)
  -s, --select      Activation (True) or deactivation (False) of select module, which filters PVC from proteome1 (default: True)
  -ss, --substitution 
                        Maximal number of compatible substitutions allowed in shared peptides alignment windows of 'minlength' size in immunity modules (default: 3)
  -tdl, --transmemb_doms_limit 
                        Parameter of select module. Proteins with trasmembrane domains >= transmemb_doms_limit are discarded (default: 3)
  -vl, --virlimit   Cut-off value for NERVirulent in the select module (0.-1) (default: 0.5)
  -vir, --virulent  Activation (True) or deactivation (False) of NERVirulent module, predictor of the probability of being a virulence factor (default: False)
  -ep, --epitopes   Activate or deactivate epitopes module (default: True)
  -m1l, --mhci_length 
                        mhci binders length (9, 10, 11 are available) (default: 9)
  -m2l, --mhcii_length 
                        mhcii binders length (9, 11, 12, 15 are available) (default: 11)
  -m1ovr, --mhci_overlap 
                        mhci-epitope overlap (default: 1)
  -m2ovr, --mhcii_overlap 
                        mhcii-epitope overlap (default: 1)
  -prt, --epitope_percentile 
                        percentile decision threshold on whick to predict epitopes from full length proteins (default: 0.9)
  -wd, --working_dir 
                        Path to working directory. If not existing, a working directory with the given path is created (default: ./)
  -nd, --NERVE_dir  Path to NERVE repository folder (download from: https://github.com/nicolagulmini/NERVE) (default: /usr/nerve_python/NERVE)
  -id, --iFeature_dir 
                        Path to iFeature repository folder (download from: https://github.com/Superzchen/iFeature) (default: /usr/nerve_python/assets/iFeature)
  -dfd, --DeepFri_dir 
                        Path to DeepFri folder (download from: https://github.com/flatironinstitute/DeepFRI) (default: /usr/nerve_python/assets/DeepFri)
  -epp, --ep_plots  Epitopes plotting script (default: True)
  ```

## Common usage and tips:

With help, you can visualize all possible arguments, which include all user-settable parameters and activation or deactivation of NERVE components.
Definitions and setting options are shown in detail for each of them here above, in the ```help``` output.

**Mandatory arguments:** 
* ```-p1  ```,```--proteome1  ``` 
* ```-g ```,  ```--gram ``` 

All the other arguments showed are optional. So, if they're not set by the user, their related default value is used during NERVE computation.

**Here is a list of different examples to show how to correctly set NERVE arguments:**

```
./NERVE.sh -p1 test.fasta -g n
```
1) This is a run example with the minimum number of arguments (the two mandatory ones). Here we're considering a FASTA file (test.fasta) with Gram negative proteins.

```
./NERVE.sh -p1 test.fasta -g n -s False -vir True -vl 0.8
```
2) The same file with Gram negative proteins will be analyzed activating NERVirulent and setting its cut-off to 0.8 (instead of default " 0.5 "). Without activating Select (```-s False```), PVCs are not filtered. So, in the results file, there will be the list of all proteins with their related analyzed info.

```
./NERVE.sh -p1 test.fasta -g n -p2 test2.fasta
```
3) In this case, adding test2.fasta, the component Conservation is automatically activated, allowing the user to infer antigen conservation. 
For more info about components working visit the FAQ section of https://nerve-bio.org

