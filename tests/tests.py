#!/usr/bin/python3
"""Tests NERVE"""

import os, shutil
from subprocess import getstatusoutput

# global test variables

PRG = "NERVE.sh"
NERVE_VERSION="v0.0.6"
RUN = f"docker run --network nerve-network -p 8880:8880 -i -v $(pwd):/workdir nerve:{NERVE_VERSION}"

WORKDIR = "./tests/output_data/"
INPUTDIR = "./tests/input_data/"
TEST1 = ('proteome1', f'--gram n --proteome1 {INPUTDIR}proteome1.fasta --working_dir\
          {WORKDIR}proteome1 --annotation False --razor False --select False --mouse False --epitopes False')
TEST2 = ('proteome1', f'--gram n --proteome1 {INPUTDIR}proteome1.fasta --working_dir\
          {WORKDIR}proteome1 --annotation True --razor False --select False --mouse False --epitopes False')
TEST3 = ('proteome1', f'--gram n --proteome1 {INPUTDIR}proteome1.fasta --working_dir\
          {WORKDIR}proteome1 --annotation True --razor False --select True --mouse False --epitopes False')
TEST4 = ('proteome1', f'--gram n --proteome1 {INPUTDIR}proteome1.fasta --working_dir\
          {WORKDIR}proteome1 --annotation True --razor False --select True --mouse True --epitopes False')
TEST5 = ('proteome2', f'--gram n --proteome1 {INPUTDIR}proteome2.fasta --working_dir\
          {WORKDIR}proteome2 --annotation True --razor True --select True --mouse True --epitopes False')

TEST6 = ('proteome3', f'--gram p --proteome1 {INPUTDIR}proteome3.fasta --working_dir\
          {WORKDIR}proteome3 --annotation True --razor True --select True --mouse True --epitopes False')
TEST7 = ('proteome4', f'--gram p --proteome1 {INPUTDIR}proteome4.fasta --working_dir\
          {WORKDIR}proteome4 --annotation True --razor True --select True --mouse True --epitopes False')
TEST8 = ('proteome4', f'--gram p --proteome1 {INPUTDIR}proteome4.fasta \
         -p2 {INPUTDIR}proteome5.fasta \
         --working_dir {WORKDIR}proteome4 --annotation True --razor True --select \
         True --mouse True --epitopes False')

TEST9 = ('proteome4', f'--gram p --proteome1 {INPUTDIR}proteome4.fasta --working_dir\
          {WORKDIR}proteome4 --annotation True --razor True --select True --mouse True --epitopes True', 
          [
              "epitope/P0C277/heatmap_pbs_MHC1_P0C277.png",
              "epitope/P0C277/MHC1_epitopes_FILTERED_P0C277.csv",
              "epitope/P0C277/Promiscuous_binders_MHC1_P0C277.csv",
              "epitope/P0C277/mhci_epitopes_P0C277.csv",
              "epitope/P0C277/mhcii_epitopes_P0C277.csv"
          ]
)

TEST10 = ('proteome6', f'--gram p --proteome1 {INPUTDIR}proteome6.fasta \
            --working_dir {WORKDIR}proteome6 --annotation True --razor True --select \
            True --mouse True --epitopes False', """>tr|Q0PC84|Q0PC84_CAMJE Basal-body rod modification protein FlgD OS=Campylobacter jejuni subsp. jejuni serotype O:2 (strain ATCC 700819 / NCTC 11168) OX=192222 GN=flgD PE=3 SV=1
MISSSDWNLNTTATTSGTTSS@GSTSGTTRTDSSSSSGIVSNPNATLDKDAFLKLLLIELQHQDPTDPMDSDKMLTQTSQLSALEMQQNTNTTMQKMVETMQKLSDSFSTSMSTSALGAIGKMATVSDNKIKLTGADELIALKMYLPEDSDENGVTLEIYDSNNKLVFSEKSDAKSISQGLFTMEWPGRNNDGVYAGDGEYTVKMVYNNKNGEKITANYGTYPIEGVVFKDGVAYAKMAGQEVPFDAIQEITDYKLGSSSSTGGSGSSGDSSGGSSDGDSSGSGSTEDGDKEEKA
""")

TEST11 = ('proteome7', f'--gram p --proteome1 {INPUTDIR}proteome7.fasta \
            --working_dir {WORKDIR}proteome7 --annotation True --razor True --select \
            True --mouse True --epitopes False', """>*tr|Q0PC84|Q0PC84_CAMJE Basal-body rod modification protein FlgD OS=Campylobacter jejuni subsp. jejuni serotype O:2 (strain ATCC 700819 / NCTC 11168) OX=192222 GN=flgD PE=3 SV=1
MISSSDWNLNTTATTSGTTSSGSTSGTTRTDSSSSSGIVSNPNATLDKDAFLKLLLIELQHQDPTDPMDSDKMLTQTSQLSALEMQQNTNTTMQKMVETMQKLSDSFSTSMSTSALGAIGKMATVSDNKIKLTGADELIALKMYLPEDSDENGVTLEIYDSNNKLVFSEKSDAKSISQGLFTMEWPGRNNDGVYAGDGEYTVKMVYNNKNGEKITANYGTYPIEGVVFKDGVAYAKMAGQEVPFDAIQEITDYKLGSSSSTGGSGSSGDSSGGSSDGDSSGSGSTEDGDKEEKAC
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

def test_arg():
    """Uses command-line arguments"""
    for proteome1, arg in [TEST1, TEST2, TEST3, TEST4, TEST5, TEST6, TEST7, TEST8]:
        rv, out = getstatusoutput(f'{RUN} {arg}')
        print("rv:\n", rv)
        print("out:\n", out)
        assert rv == 0
        assert out.endswith("End NERVE computation successfully.")
        assert os.path.isfile(os.path.join(WORKDIR, proteome1, "vaccine_candidates.csv"))
        assert os.path.isfile(os.path.join(WORKDIR, proteome1, "discarded_proteins.csv"))

def test_files():
    """Tests epitope output files"""
    for proteome1, arg, expected1 in [TEST9]:
        rv, out = getstatusoutput(f'{RUN} {arg}')
        print("rv:\n", rv)
        print("out:\n", out)
        for file in expected1:
            assert os.path.isfile(os.path.join(WORKDIR, proteome1, file)) == True

def test_exception_1():
    """Test presence of input sequences with wrong format"""
    for proteome1, arg, expected1 in [TEST10]:
        rv, out = getstatusoutput(f'{RUN} {arg}')
        print("rv:\n", rv)
        print("out:\n", out)
        print("expected1:\n", expected1, "\neffective1:\n", open(os.path.join(WORKDIR, proteome1, 
                                    f"discarded_sequences_{proteome1}.fasta"), 'r').read())
        assert expected1 == open(os.path.join(WORKDIR, proteome1, 
                                              f"discarded_sequences_{proteome1}.fasta"), 'r').read()
    
def test_exception_2():
    """Test presence of input sequences with wrong format"""
    for proteome1, arg, expected1 in [TEST11]:
        rv, out = getstatusoutput(f'{RUN} {arg}')
        print("rv:\n", rv)
        print("out:\n", out)
        print("expected1:\n", expected1, "\neffective1:\n", open(os.path.join(WORKDIR, proteome1, 
                                    f"cleaned_{proteome1}.fasta"), 'r').read())
        assert expected1 == open(os.path.join(WORKDIR, proteome1, 
                                              f"cleaned_{proteome1}.fasta"), 'r').read()