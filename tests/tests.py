#!python3

"""Tests NERVE"""

import os, shutil
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

TEST6 = ('proteome3', f'--gram p --proteome1 {INPUTDIR}proteome3.fasta --working_dir\
          {WORKDIR}proteome3 --annotation True --razor True --select True --mouse True', """,id ,uniprot_accession_code,length,transmembrane_doms,localization,localization score,virulence_probability,adhesin_probability,conservation_score,list_of_shared_human_peps,list_of_shared_mouse_peps,list_of_shared_conserv_proteome_peps,human_peptides_sum,mouse_peptides_sum,annotations,list_of_peptides_from_comparison_with_mhcpep_sapiens,list_of_peptides_from_comparison_with_mhcpep_mouse,sequence,original_sequence_if_razor,tmhmm_seq
""", """,id ,uniprot_accession_code,length,transmembrane_doms,localization,localization score,virulence_probability,adhesin_probability,conservation_score,list_of_shared_human_peps,list_of_shared_mouse_peps,list_of_shared_conserv_proteome_peps,human_peptides_sum,mouse_peptides_sum,annotations,list_of_peptides_from_comparison_with_mhcpep_sapiens,list_of_peptides_from_comparison_with_mhcpep_mouse,sequence,original_sequence_if_razor,tmhmm_seq
0,tr|N9Y0T0|N9Y0T0_9CLOT Coenzyme A biosynthesis bifunctional protein CoaBC OS=Clostridium thermobutyricum OX=29372 GN=coaBC PE=3 SV=1,N9Y0T0,394,0,Cytoplasmic,7.5,0.1401,0.0358,,19,19,0,0.1701,0.1497,DeepFri predictions: ,,,MYEDKCVVIGVTGGIAVYKALDVISALRKKGVKTKVIMTESATKFVTPLTFQSISQNMVITDMFAEPKAWKIQHISLAQEADIMLIAPATANVIGKVANGIADDMLTTTIMATKAKVIFSPAMNTNMYENPIVQENIEKLKRFGYEFIEPDSGRLACGDIGKGKLPKPEVIVDEVLKNLYPKKDLVGKKVLVTAGPTKAPIDPVRYITNRSTGKMGYAIATEARDRGAEVTLVSGADFLECPKGIELINVQTNSEMREEVLKHYKNSDIVIKSAAVADYKPKNYSKEKIKKTGDDLKLELERDNDILLELGKLKKNQILVGFAAESNSLLENANKKLKNKNLDFIVANDITSKETGFGSDNNKVFIISKNGKVDELDTMSKRDVARSIFDNILE,,
""")
TEST7 = ('proteome4', f'--gram p --proteome1 {INPUTDIR}proteome4.fasta --working_dir\
          {WORKDIR}proteome4 --annotation True --razor True --select True --mouse True', """,id ,uniprot_accession_code,length,transmembrane_doms,localization,localization score,virulence_probability,adhesin_probability,conservation_score,list_of_shared_human_peps,list_of_shared_mouse_peps,list_of_shared_conserv_proteome_peps,human_peptides_sum,mouse_peptides_sum,annotations,list_of_peptides_from_comparison_with_mhcpep_sapiens,list_of_peptides_from_comparison_with_mhcpep_mouse,sequence,original_sequence_if_razor,tmhmm_seq
0,tr|Q0PA76|Q0PA76_CAMJE Uncharacterized protein OS=Campylobacter jejuni subsp. jejuni serotype O:2 (strain ATCC 700819 / NCTC 11168) OX=192222 GN=Cj0814 PE=4 SV=1,Q0PA76,251,0,Unknown,0.0,0.9651,0.992,,0,0,0,0.0,0.0,DeepFri predictions: ,,,MITQTMQSKESKESKENSKISFANAFLKQNASKLNEIQNANSQTLARSEALNSTNTTNTSNNTNFSISSKTSSPNYDISSEFKNSIYTLKYKQVDISNTSTNTAYGYSVDKDGYMGSDFNKAAGLPEDFKIHKSTLDEIKKAAENDPVVSSTKEYLGVSSYYSNIDIANTIKQYYNLFSNALGQSFSNDKTSFSEADINSMPSGYGVSGTQWMDFNEPSNRMNITGLKDFSNSLISNVYKTPEQAKEADEI,,
1,tr|Q0PC84|Q0PC84_CAMJE Basal-body rod modification protein FlgD OS=Campylobacter jejuni subsp. jejuni serotype O:2 (strain ATCC 700819 / NCTC 11168) OX=192222 GN=flgD PE=3 SV=1,Q0PC84,294,0,Unknown,0.0,0.9646,0.9789,,0,0,0,0.0,0.0,DeepFri predictions: ,,,MISSSDWNLNTTATTSGTTSSGSTSGTTRTDSSSSSGIVSNPNATLDKDAFLKLLLIELQHQDPTDPMDSDKMLTQTSQLSALEMQQNTNTTMQKMVETMQKLSDSFSTSMSTSALGAIGKMATVSDNKIKLTGADELIALKMYLPEDSDENGVTLEIYDSNNKLVFSEKSDAKSISQGLFTMEWPGRNNDGVYAGDGEYTVKMVYNNKNGEKITANYGTYPIEGVVFKDGVAYAKMAGQEVPFDAIQEITDYKLGSSSSTGGSGSSGDSSGGSSDGDSSGSGSTEDGDKEEKA,,
""", """,id ,uniprot_accession_code,length,transmembrane_doms,localization,localization score,virulence_probability,adhesin_probability,conservation_score,list_of_shared_human_peps,list_of_shared_mouse_peps,list_of_shared_conserv_proteome_peps,human_peptides_sum,mouse_peptides_sum,annotations,list_of_peptides_from_comparison_with_mhcpep_sapiens,list_of_peptides_from_comparison_with_mhcpep_mouse,sequence,original_sequence_if_razor,tmhmm_seq
0,"tr|Q0PAN8|Q0PAN8_CAMJE Putative DNA polymerase III, delta subunit OS=Campylobacter jejuni subsp. jejuni serotype O:2 (strain ATCC 700819 / NCTC 11168) OX=192222 GN=holA PE=4 SV=1",Q0PAN8,321,0,Cytoplasmic,7.5,0.8001,0.2387,,0,0,0,0.0,0.0,DeepFri predictions: ,,,MYRKELQILLSKDSIPNFFFLYGADNFQSELYAEFIKEKYKPDETLKLFFEEYNFTRASDFLSAGSLFSEKKLLEIKTSKKIPTKDLKVLVELCKNNTDNFFLLELYDESSKQSDIEKIFSPHFVRFFKANGAKEGVELLSIKAKQLGVEITQNALFTLFTSFDENLYLAASELNKFSGLRVDEKTIEQYCYSLNTGSFESFFDKILKKQDFKSELEKILDNFNEIALINSLYNSFYRLFKIALYAKINGKIDFKELLGYTPPPQVGQNLSSQAFSLKIEQYKEIFTLLLKSEYELKTNPKLVKKEFLISNLLKLARILKN,,
1,tr|Q0PAP1|Q0PAP1_CAMJE Hydrogenase isoenzymes formation protein OS=Campylobacter jejuni subsp. jejuni serotype O:2 (strain ATCC 700819 / NCTC 11168) OX=192222 GN=hypE PE=3 SV=1,Q0PAP1,324,0,Unknown,0.0,0.1983,0.0169,,0,0,0,0.0,0.0,DeepFri predictions: ,,,MKNISLAHGGGGEEMNELLTKLFKIFDNEILNANNDAAILGNLALSTDSFVLSPIFLDEEVNIGKLCVCGSINDVLMVGAKPKYLSLGLILEEGFELEKLERILKSIKEECEKCGVMLVCGDTKVVPKGKADEIYINTTALGEIISKKESKNIKAGLSILLSGDIGRHGASVLIKRNELEADIKSDCKALNKEVLELLEKDIKVVAMRDATRGGLSAVLNEWAKQSGNDLLIFEEKIIVQDEVLGLCELFGYEAYELANEGTFILCVEKEDELKALEILKKYNVNASIIGEVLEEKKARVILQNAYGAKRFLESPKGELLPRIC,,
""")
TEST8 = ('proteome4', f'--gram p --proteome1 {INPUTDIR}proteome4.fasta \
         -p2 {INPUTDIR}proteome5.fasta \
         --working_dir {WORKDIR}proteome4 --annotation True --razor True --select \
         True --mouse True ', """,id ,uniprot_accession_code,length,transmembrane_doms,localization,localization score,virulence_probability,adhesin_probability,conservation_score,list_of_shared_human_peps,list_of_shared_mouse_peps,list_of_shared_conserv_proteome_peps,human_peptides_sum,mouse_peptides_sum,annotations,list_of_peptides_from_comparison_with_mhcpep_sapiens,list_of_peptides_from_comparison_with_mhcpep_mouse,sequence,original_sequence_if_razor,tmhmm_seq
0,tr|Q0PA76|Q0PA76_CAMJE Uncharacterized protein OS=Campylobacter jejuni subsp. jejuni serotype O:2 (strain ATCC 700819 / NCTC 11168) OX=192222 GN=Cj0814 PE=4 SV=1,Q0PA76,251,0,Unknown,0.0,0.9651,0.992,0.0,0,0,0,0.0,0.0,DeepFri predictions: ,,,MITQTMQSKESKESKENSKISFANAFLKQNASKLNEIQNANSQTLARSEALNSTNTTNTSNNTNFSISSKTSSPNYDISSEFKNSIYTLKYKQVDISNTSTNTAYGYSVDKDGYMGSDFNKAAGLPEDFKIHKSTLDEIKKAAENDPVVSSTKEYLGVSSYYSNIDIANTIKQYYNLFSNALGQSFSNDKTSFSEADINSMPSGYGVSGTQWMDFNEPSNRMNITGLKDFSNSLISNVYKTPEQAKEADEI,,
1,tr|Q0PC84|Q0PC84_CAMJE Basal-body rod modification protein FlgD OS=Campylobacter jejuni subsp. jejuni serotype O:2 (strain ATCC 700819 / NCTC 11168) OX=192222 GN=flgD PE=3 SV=1,Q0PC84,294,0,Unknown,0.0,0.9646,0.9789,1.0,0,0,286,0.0,0.0,DeepFri predictions: ,,,MISSSDWNLNTTATTSGTTSSGSTSGTTRTDSSSSSGIVSNPNATLDKDAFLKLLLIELQHQDPTDPMDSDKMLTQTSQLSALEMQQNTNTTMQKMVETMQKLSDSFSTSMSTSALGAIGKMATVSDNKIKLTGADELIALKMYLPEDSDENGVTLEIYDSNNKLVFSEKSDAKSISQGLFTMEWPGRNNDGVYAGDGEYTVKMVYNNKNGEKITANYGTYPIEGVVFKDGVAYAKMAGQEVPFDAIQEITDYKLGSSSSTGGSGSSGDSSGGSSDGDSSGSGSTEDGDKEEKA,,
""", """,id ,uniprot_accession_code,length,transmembrane_doms,localization,localization score,virulence_probability,adhesin_probability,conservation_score,list_of_shared_human_peps,list_of_shared_mouse_peps,list_of_shared_conserv_proteome_peps,human_peptides_sum,mouse_peptides_sum,annotations,list_of_peptides_from_comparison_with_mhcpep_sapiens,list_of_peptides_from_comparison_with_mhcpep_mouse,sequence,original_sequence_if_razor,tmhmm_seq
0,"tr|Q0PAN8|Q0PAN8_CAMJE Putative DNA polymerase III, delta subunit OS=Campylobacter jejuni subsp. jejuni serotype O:2 (strain ATCC 700819 / NCTC 11168) OX=192222 GN=holA PE=4 SV=1",Q0PAN8,321,0,Cytoplasmic,7.5,0.8001,0.2387,0.0,0,0,0,0.0,0.0,DeepFri predictions: ,,,MYRKELQILLSKDSIPNFFFLYGADNFQSELYAEFIKEKYKPDETLKLFFEEYNFTRASDFLSAGSLFSEKKLLEIKTSKKIPTKDLKVLVELCKNNTDNFFLLELYDESSKQSDIEKIFSPHFVRFFKANGAKEGVELLSIKAKQLGVEITQNALFTLFTSFDENLYLAASELNKFSGLRVDEKTIEQYCYSLNTGSFESFFDKILKKQDFKSELEKILDNFNEIALINSLYNSFYRLFKIALYAKINGKIDFKELLGYTPPPQVGQNLSSQAFSLKIEQYKEIFTLLLKSEYELKTNPKLVKKEFLISNLLKLARILKN,,
1,tr|Q0PAP1|Q0PAP1_CAMJE Hydrogenase isoenzymes formation protein OS=Campylobacter jejuni subsp. jejuni serotype O:2 (strain ATCC 700819 / NCTC 11168) OX=192222 GN=hypE PE=3 SV=1,Q0PAP1,324,0,Unknown,0.0,0.1983,0.0169,0.0,0,0,0,0.0,0.0,DeepFri predictions: ,,,MKNISLAHGGGGEEMNELLTKLFKIFDNEILNANNDAAILGNLALSTDSFVLSPIFLDEEVNIGKLCVCGSINDVLMVGAKPKYLSLGLILEEGFELEKLERILKSIKEECEKCGVMLVCGDTKVVPKGKADEIYINTTALGEIISKKESKNIKAGLSILLSGDIGRHGASVLIKRNELEADIKSDCKALNKEVLELLEKDIKVVAMRDATRGGLSAVLNEWAKQSGNDLLIFEEKIIVQDEVLGLCELFGYEAYELANEGTFILCVEKEDELKALEILKKYNVNASIIGEVLEEKKARVILQNAYGAKRFLESPKGELLPRIC,,
""")

TEST9 = ('proteome4', f'--gram p --proteome1 {INPUTDIR}proteome4.fasta --working_dir\
          {WORKDIR}proteome4 --annotation True --razor True --select True --mouse True', ['mouse_immunity_raw_output.txt',
 'cleaned_proteome4.fasta',
 'logfile.log',
 '_MF_pred_scores.json',
 'vaccine_candidates.csv',
 'payload.json',
 'discarded_sequences_proteome4.fasta',
 'autoimmunity_raw_output.txt',
 '_MF_predictions.csv',
 'discarded_proteins.csv']
)

TEST10 = ('proteome6', f'--gram p --proteome1 {INPUTDIR}proteome6.fasta \
            --working_dir {WORKDIR}proteome6 --annotation True --razor True --select \
            True --mouse True', """>tr|Q0PC84|Q0PC84_CAMJE Basal-body rod modification protein FlgD OS=Campylobacter jejuni subsp. jejuni serotype O:2 (strain ATCC 700819 / NCTC 11168) OX=192222 GN=flgD PE=3 SV=1
MISSSDWNLNTTATTSGTTSSGSTSGTTRTDSSSSSGIVSNPNATLDKDAFLKLLLIELQHQDPTDPMDSDKMLTQTSQLSALEMQQNTNTTMQKMVETMQKLSDSFSTSMSTSALGAIGKMATVSDNKIKLTGADELIALKMYLPEDSDENGVTLEIYDSNNKLVFSEKSDAKSISQGLFTMEWPGRNNDGVYAGDGEYTVKMVYNNKNGEKITANYGTYPIEGVVFKDGVAYAKMAGQEVPFDAIQEITDYKLGSSSSTGGSGSSGDSSGGSSDGDSSGSGSTEDGDKEEKA
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
    for proteome1, arg, expected1, expected2 in [TEST1, TEST2, TEST3, TEST4, TEST5, TEST6, TEST7, TEST8]:
        rv, out = getstatusoutput(f'{RUN} {arg}')
        print("rv:\n", rv)
        print("out:\n", out)
        print("expected1:\n", expected1, "effective1:\n", open(os.path.join(WORKDIR, proteome1, "vaccine_candidates.csv"), 'r').read())
        print("expected2:\n", expected2, "effective2:\n", open(os.path.join(WORKDIR, proteome1, "discarded_proteins.csv"), 'r').read())
        assert rv == 0
        assert out.endswith("End NERVE computation successfully.")
        assert open(os.path.join(WORKDIR, proteome1, "vaccine_candidates.csv"), 'r').read() == expected1
        assert open(os.path.join(WORKDIR, proteome1, "discarded_proteins.csv"), 'r').read() == expected2
        shutil.rmtree(os.path.join(WORKDIR, proteome1))

def test_files():
    """Tests output files"""
    for proteome1, arg, expected1 in [TEST9]:
        rv, out = getstatusoutput(f'{RUN} {arg}')
        print("rv:\n", rv)
        print("out:\n", out)
        print("expected1:\n", expected1, "\neffective1:\n", os.listdir(os.path.join(WORKDIR, proteome1)))
        assert expected1 == os.listdir(os.path.join(WORKDIR, proteome1))
        shutil.rmtree(os.path.join(WORKDIR, proteome1))
        
def test_exceptions():
    """Test presence of input sequences with wrong format"""
    for proteome1, arg, expected1 in [TEST10]:
        rv, out = getstatusoutput(f'{RUN} {arg}')
        print("rv:\n", rv)
        print("out:\n", out)
        print("expected1:\n", expected1, "\neffective1:\n", open(os.path.join(WORKDIR, proteome1, "discarded_sequences_proteome6.fasta"), 'r').read())
        assert expected1 == open(os.path.join(WORKDIR, proteome1, "discarded_sequences_proteome6.fasta"), 'r').read()
        shutil.rmtree(os.path.join(WORKDIR, proteome1)) 

def test_errors():
    """Test wrong input errors"""
    print()
    
    