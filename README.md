# NERVE

Software for Reverse Vaccinology. For web-tool please visit: https://nerve-bio.org

### Stand-alone usage:
NERVE can be used as a stand alone verison taking advantage of [Docker](https://www.docker.com/) in linux systems.

1) install Docker following [these instructions](https://docs.docker.com/engine/install/) and [the post-installation procedure](https://docs.docker.com/engine/install/linux-postinstall/)
2) clone the repository:
```
git clone git@github.com:nerve-bio/NERVE.git
```
4) navigate to the correct folder:
```
cd NERVE
```
4) Run NERVE (the first time it will take a few minutes)
```
./NERVE.sh --help
```
Expected output:
```
usage: NERVE.py [-h] [-a] [-ev] -g [-ml] [-mm] [-m] [-mpsl] -p1 [-p2] [-paefilter] [-pacfilter] [-pl] [-rz] [-rl] [-s] [-ss] [-tdl] [-vl] [-vir]
                [-ep] [-m1l] [-m2l] [-m1ovr] [-m2ovr] [-prt] [-wd] [-nd] [-id] [-dfd] [-epp]

Run vaccine candidate prediction

optional arguments:
  -h, --help            show this help message and exit
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
  -paefilter, --p_ad_extracellular_filter 
                        Parameter of select module. Extracellular proteins with a probability of adhesin (pad) lower than p_ad_extracellular_filter are discarded (0.-1) (default:
                        0.38)
  -pacfilter, --p_ad_no_citoplasm_filter 
                        Parameter of select module. Non-cytoplasmic Proteins with a probability of adhesin (pad) lower than p_ad_no_citoplasm_filter are discarded (0.-1) (default:
                        0.46)
  -pl, --padlimit   Set the probability of adhesin (pad) value cut-off for proteins with 'Unknown' localization in the select module. Thus, these proteins with a pad value < cut-
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

