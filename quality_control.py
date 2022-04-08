from Bio import SeqIO
from Bio.Seq import Seq

def quality_control(path_to_fasta:str, verbosity:int)->None:
    """
    Remove sequences with non-canonical aminoacid symbols. U (Se-Cys) is substituted with C (Cys). Returns
    "non_filtered_"+input_filename and overwrites input.
    param: path_to_fasta: full path to fasta file containing the proteome without .fasta extension;
    param: verbosity: verbosity level. Print status messages if > 0:
    """
    aa_dic = {'C': 'C', 'D': 'D', 'S': 'S', 'Q': 'Q', 'K': 'K', 'I': 'I', 'P': 'P', 'T': 'T', 'F': 'F', 'N': 'N', 
              'G': 'G', 'H': 'H', 'L': 'L', 'R': 'R', 'W': 'W', 'A': 'A', 'V': 'V', 'E': 'E', 'Y': 'Y', 'M': 'M', 
              'U':'C'}
    outlist = []
    #path_to_fasta = ".".join([path_to_fasta, ".fasta"])
    fasta_list = list(SeqIO.parse(path_to_fasta, "fasta")) 
    
    filename = open("non_filtered_"+path_to_fasta.split('/')[-1], 'w')
    SeqIO.write(fasta_list, filename, "fasta")
    filename.close()
    
    for record in fasta_list:
        flag = True
        new_seq =''
        for aa in str(record.seq):
            if aa not in aa_dic:
                flag = False
                if verbosity > 0:
                    print(f'Found non-canonical aminoacid named {aa} in sequence {record.id}')
            else:
                new_seq += aa_dic[aa]
        record.seq = Seq(new_seq)
        if flag == True:
            outlist.append(record)
        elif verbosity > 0:
            print(f'Sequence {record.id} has been discarded for the presence of non-canonical aminoacids.')
        
    # filtered
    filename = open(path_to_fasta, 'w')
    SeqIO.write(outlist, filename, "fasta")
    filename.close()
        
    return None
        
