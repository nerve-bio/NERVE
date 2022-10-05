#!python3

"""Tests NERVE"""

import os
from subprocess import getstatusoutput

# global test variables

PRG = "/NERVE/NERVE/NERVE_main.py"
RUN = f'python {PRG}'
WORKDIR = "/my_data/tests/output_data/"
INPUTDIR = "/my_data/tests/input_data/"
TEST1 = ('proteome1', f'--gram n --proteome1 {INPUTDIR}proteome1.fasta --working_dir\
          {WORKDIR}proteome1 --annotation False --razor False --select False --mouse False', """,id ,uniprot_accession_code,length,transmembrane_doms,localization,localization score,virulence_probability,adhesin_probability,conservation_score,list_of_shared_human_peps,list_of_shared_mouse_peps,list_of_shared_conserv_proteome_peps,human_peptides_sum,mouse_peptides_sum,annotations,list_of_peptides_from_comparison_with_mhcpep_sapiens,list_of_peptides_from_comparison_with_mhcpep_mouse,sequence,original_sequence_if_razor,tmhmm_seq
0,sp|P0C277|FOLD_NEIMB Bifunctional protein FolD OS=Neisseria meningitidis serogroup B (strain MC58) OX=122586 GN=folD PE=3 SV=1,P0C277,284,0,Periplasmic,5.96,0.1789,0.0243,,117,0,0,1.0599,0,,,,MSAQLINGKEVSQKRLQAVAEAVAQRQQNNLHHPCLAVVLVGGDPASAVYVRNKKTACQKCGIKSLSYELPESTSQEELLALVDRLNADSEVDGILVQLPLPKHLDSQAVLERISPDKDVDGFHPYNVGRLAVKMPLMRPCTPKGVMTLLEAYGIDPKGKKAVVVGASNIVGRPQALELLLARATVTVCHSATENLTDEVAGADILVVGVGIPNFVKGEWIKPGAVVIDVGINRLDDGSLCGDVEFETAKERAAMITPVPGGVGPMTIATLMENTLHAASLHDA,,
""", """,id ,uniprot_accession_code,length,transmembrane_doms,localization,localization score,virulence_probability,adhesin_probability,conservation_score,list_of_shared_human_peps,list_of_shared_mouse_peps,list_of_shared_conserv_proteome_peps,human_peptides_sum,mouse_peptides_sum,annotations,list_of_peptides_from_comparison_with_mhcpep_sapiens,list_of_peptides_from_comparison_with_mhcpep_mouse,sequence,original_sequence_if_razor,tmhmm_seq
""")
TEST2 = ('proteome1', f'--gram n --proteome1 {INPUTDIR}proteome1.fasta --working_dir\
          {WORKDIR}proteome1 --annotation True --razor False --select False --mouse False', """,id ,uniprot_accession_code,length,transmembrane_doms,localization,localization score,virulence_probability,adhesin_probability,conservation_score,list_of_shared_human_peps,list_of_shared_mouse_peps,list_of_shared_conserv_proteome_peps,human_peptides_sum,mouse_peptides_sum,annotations,list_of_peptides_from_comparison_with_mhcpep_sapiens,list_of_peptides_from_comparison_with_mhcpep_mouse,sequence,original_sequence_if_razor,tmhmm_seq
0,sp|P0C277|FOLD_NEIMB Bifunctional protein FolD OS=Neisseria meningitidis serogroup B (strain MC58) OX=122586 GN=folD PE=3 SV=1,P0C277,284,0,Periplasmic,5.96,0.1789,0.0243,,117,0,0,1.0599,0,"DeepFri predictions: cyclohydrolase activity | oxidoreductase activity, acting on the CH-NH group of donors, NAD or NADP as acceptor | oxidoreductase activity, acting on the CH-NH group of donors | hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds, in cyclic amidines | hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds",,,MSAQLINGKEVSQKRLQAVAEAVAQRQQNNLHHPCLAVVLVGGDPASAVYVRNKKTACQKCGIKSLSYELPESTSQEELLALVDRLNADSEVDGILVQLPLPKHLDSQAVLERISPDKDVDGFHPYNVGRLAVKMPLMRPCTPKGVMTLLEAYGIDPKGKKAVVVGASNIVGRPQALELLLARATVTVCHSATENLTDEVAGADILVVGVGIPNFVKGEWIKPGAVVIDVGINRLDDGSLCGDVEFETAKERAAMITPVPGGVGPMTIATLMENTLHAASLHDA,,
""", """,id ,uniprot_accession_code,length,transmembrane_doms,localization,localization score,virulence_probability,adhesin_probability,conservation_score,list_of_shared_human_peps,list_of_shared_mouse_peps,list_of_shared_conserv_proteome_peps,human_peptides_sum,mouse_peptides_sum,annotations,list_of_peptides_from_comparison_with_mhcpep_sapiens,list_of_peptides_from_comparison_with_mhcpep_mouse,sequence,original_sequence_if_razor,tmhmm_seq
""")
TEST3 = ('proteome1', f'--gram n --proteome1 {INPUTDIR}proteome1.fasta --working_dir\
          {WORKDIR}proteome1 --annotation True --razor False --select True --mouse False', """,id ,uniprot_accession_code,length,transmembrane_doms,localization,localization score,virulence_probability,adhesin_probability,conservation_score,list_of_shared_human_peps,list_of_shared_mouse_peps,list_of_shared_conserv_proteome_peps,human_peptides_sum,mouse_peptides_sum,annotations,list_of_peptides_from_comparison_with_mhcpep_sapiens,list_of_peptides_from_comparison_with_mhcpep_mouse,sequence,original_sequence_if_razor,tmhmm_seq
""", """,id ,uniprot_accession_code,length,transmembrane_doms,localization,localization score,virulence_probability,adhesin_probability,conservation_score,list_of_shared_human_peps,list_of_shared_mouse_peps,list_of_shared_conserv_proteome_peps,human_peptides_sum,mouse_peptides_sum,annotations,list_of_peptides_from_comparison_with_mhcpep_sapiens,list_of_peptides_from_comparison_with_mhcpep_mouse,sequence,original_sequence_if_razor,tmhmm_seq
0,sp|P0C277|FOLD_NEIMB Bifunctional protein FolD OS=Neisseria meningitidis serogroup B (strain MC58) OX=122586 GN=folD PE=3 SV=1,P0C277,284,0,Periplasmic,5.96,0.1789,0.0243,,117,0,0,1.0599,0,"DeepFri predictions: cyclohydrolase activity | oxidoreductase activity, acting on the CH-NH group of donors, NAD or NADP as acceptor | oxidoreductase activity, acting on the CH-NH group of donors | hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds, in cyclic amidines | hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds",,,MSAQLINGKEVSQKRLQAVAEAVAQRQQNNLHHPCLAVVLVGGDPASAVYVRNKKTACQKCGIKSLSYELPESTSQEELLALVDRLNADSEVDGILVQLPLPKHLDSQAVLERISPDKDVDGFHPYNVGRLAVKMPLMRPCTPKGVMTLLEAYGIDPKGKKAVVVGASNIVGRPQALELLLARATVTVCHSATENLTDEVAGADILVVGVGIPNFVKGEWIKPGAVVIDVGINRLDDGSLCGDVEFETAKERAAMITPVPGGVGPMTIATLMENTLHAASLHDA,,
""")
TEST4 = ('proteome1', f'--gram n --proteome1 {INPUTDIR}proteome1.fasta --working_dir\
          {WORKDIR}proteome1 --annotation True --razor False --select True --mouse True', """,id ,uniprot_accession_code,length,transmembrane_doms,localization,localization score,virulence_probability,adhesin_probability,conservation_score,list_of_shared_human_peps,list_of_shared_mouse_peps,list_of_shared_conserv_proteome_peps,human_peptides_sum,mouse_peptides_sum,annotations,list_of_peptides_from_comparison_with_mhcpep_sapiens,list_of_peptides_from_comparison_with_mhcpep_mouse,sequence,original_sequence_if_razor,tmhmm_seq
""", """,id ,uniprot_accession_code,length,transmembrane_doms,localization,localization score,virulence_probability,adhesin_probability,conservation_score,list_of_shared_human_peps,list_of_shared_mouse_peps,list_of_shared_conserv_proteome_peps,human_peptides_sum,mouse_peptides_sum,annotations,list_of_peptides_from_comparison_with_mhcpep_sapiens,list_of_peptides_from_comparison_with_mhcpep_mouse,sequence,original_sequence_if_razor,tmhmm_seq
0,sp|P0C277|FOLD_NEIMB Bifunctional protein FolD OS=Neisseria meningitidis serogroup B (strain MC58) OX=122586 GN=folD PE=3 SV=1,P0C277,284,0,Periplasmic,5.96,0.1789,0.0243,,117,117,0,1.0599,1.0317,"DeepFri predictions: cyclohydrolase activity | oxidoreductase activity, acting on the CH-NH group of donors, NAD or NADP as acceptor | oxidoreductase activity, acting on the CH-NH group of donors | hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds, in cyclic amidines | hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds",,,MSAQLINGKEVSQKRLQAVAEAVAQRQQNNLHHPCLAVVLVGGDPASAVYVRNKKTACQKCGIKSLSYELPESTSQEELLALVDRLNADSEVDGILVQLPLPKHLDSQAVLERISPDKDVDGFHPYNVGRLAVKMPLMRPCTPKGVMTLLEAYGIDPKGKKAVVVGASNIVGRPQALELLLARATVTVCHSATENLTDEVAGADILVVGVGIPNFVKGEWIKPGAVVIDVGINRLDDGSLCGDVEFETAKERAAMITPVPGGVGPMTIATLMENTLHAASLHDA,,
""")
TEST5 = ('proteome2', f'--gram n --proteome1 {INPUTDIR}proteome2.fasta --working_dir\
          {WORKDIR}proteome2 --annotation True --razor True --select True --mouse True', """,id ,uniprot_accession_code,length,transmembrane_doms,localization,localization score,virulence_probability,adhesin_probability,conservation_score,list_of_shared_human_peps,list_of_shared_mouse_peps,list_of_shared_conserv_proteome_peps,human_peptides_sum,mouse_peptides_sum,annotations,list_of_peptides_from_comparison_with_mhcpep_sapiens,list_of_peptides_from_comparison_with_mhcpep_mouse,sequence,original_sequence_if_razor,tmhmm_seq
0,tr|Q0P9Y0|Q0P9Y0_CAMJE Putative ABC-type amino-acid transporter permease protein OS=Campylobacter jejuni subsp. jejuni serotype O:2 (strain ATCC 700819 / NCTC 11168) OX=192222 GN=Cj0919c PE=3 SV=1,Q0P9Y0,219,3,CytoplasmicMembrane,10.0,0.6636,0.7809,,0,0,0,0.0,0.0,DeepFri predictions: ,,,GTIGFSLYTSSVMAEIIRGGLNSIPKGQFEAAYSQGFGKFFTLFYIILPQTFRKIIPALLSQIVTTVKDTAYLAGLGIAELTYNSKTILAKLTSFEEI,MENVFNAQNIEFLMQGLFLTLKIALATCIISIVFGTFLAITKNYGDRLSKFLAACYIDIFRNTPLLLWMLAACFVLPVFFGQFPQAFWGTIGFSLYTSSVMAEIIRGGLNSIPKGQFEAAYSQGFGKFFTLFYIILPQTFRKIIPALLSQIVTTVKDTAYLAGLGIAELTYNSKTILAKLTSFEEILAMIGVVAGIYFIICFSLSMLVRYYAKKTAYIS,ooooooooooooooooooooooMMMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMMMooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooMMMMMMMMMMMMMMMMMMMMMMMiiiiiiiiii
""", """,id ,uniprot_accession_code,length,transmembrane_doms,localization,localization score,virulence_probability,adhesin_probability,conservation_score,list_of_shared_human_peps,list_of_shared_mouse_peps,list_of_shared_conserv_proteome_peps,human_peptides_sum,mouse_peptides_sum,annotations,list_of_peptides_from_comparison_with_mhcpep_sapiens,list_of_peptides_from_comparison_with_mhcpep_mouse,sequence,original_sequence_if_razor,tmhmm_seq
""")

# test functions

def test_exists() -> None:
    """Tests if directory is correct"""
    assert os.path.exists(PRG)

def test_usage() -> None:
    """Prints usage"""
    for arg in ['-h', '--help']:
        rv, out = getstatusoutput(f'{RUN} {arg}')
        assert rv == 0 
        assert out.lower().startswith('usage:')

def test_arg():
    """Uses command-line arguments"""
    for proteome1, arg, expected1, expected2 in [TEST1, TEST2, TEST3, TEST4, TEST5]:
        rv, out = getstatusoutput(f'{RUN} {arg}')
        print("expected1:\n", expected1, "effective1:\n", open(os.path.join(WORKDIR, proteome1, "vaccine_candidates.csv"), 'r').read())
        print("expected2:\n", expected2, "effective2:\n", open(os.path.join(WORKDIR, proteome1, "discarded_proteins.csv"), 'r').read())
        assert rv == 0
        assert out.endswith("End NERVE computation successfully.")
        assert open(os.path.join(WORKDIR, proteome1, "vaccine_candidates.csv"), 'r').read() == expected1
        assert open(os.path.join(WORKDIR, proteome1, "discarded_proteins.csv"), 'r').read() == expected2