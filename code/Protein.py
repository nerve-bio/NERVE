#!/usr/local/bin/python
"""Class and methods to store entry information"""

import numpy as np

class Protein:
		
	def __init__(self, identifier, sequence_string):
		self.id = identifier
		self.score = None        
		self.accession = identifier.split('|')[1] if identifier.count("|")==2 else None
		self.sequence = sequence_string # the sequence used for the analyses
		self.original_sequence_if_razor = None # put the original sequence if razor is performed
		self.sequence_out = None
		self.length = len(sequence_string)
		self.localization = None
		self.p_ad = None
		self.transmembrane_doms = None
		self.tmhmm_seq = None
		self.list_of_shared_human_peps = [] # in a dictionary with match, query and starting position (start_pos)
		self.list_of_shared_mouse_peps = []
		self.list_of_shared_conserv_proteome_peps = []
		self.list_of_peptides_from_comparison_with_mhcpep_sapiens = [] # here a list of mhcpep match
		self.list_of_peptides_from_comparison_with_mhcpep_mouse = [] # here a list of mhcpep match
		self.razor_loops = []
		self.p_vir = None
		self.sapiens_peptides_sum = None
		self.mouse_peptides_sum = None
		self.conservation_score = None
		self.annotations = None
		self.model_raw_data = []
		self.model_raw_data1 = []
		self.MHC1_binders = []
		self.MHC2_binders = []
		self.MHC1_pb_binders = []
		self.MHC2_pb_binders = []
		

	def print_information(self):
		print("Information about protein " + str(self.id) + ":")
		print('	accession:', self.accession)
		print('	length of the sequence:', self.length)
		print('	localization:', self.localization)
		print('	P_ad:', self.p_ad)
		print('	transmembrane doms:', self.transmembrane_doms)
		print('	number of shared human peps:', len(self.list_of_shared_human_peps))
		print('	number of shared mouse peps:', len(self.list_of_shared_mouse_peps))
		print('	number of conservation peps:', len(self.list_of_shared_conserv_proteome_peps))
		print('	number of peps from comparison with mhcpep sapiens:', len(self.list_of_peptides_from_comparison_with_mhcpep_sapiens))
		print('	number of peps from comparison with mhcpep mouse:', len(self.list_of_peptides_from_comparison_with_mhcpep_mouse))
		print('	razor loops:', self.razor_loops)
		print('	P_vir:', self.p_vir)
		print('	autoimmunity score:', self.sapiens_peptides_sum)
		print('	mouse immunity score:', self.mouse_peptides_sum)
		print('	conservation score:', self.conservation_score)
		print('	annotations:', self.annotations)
		      
	def standardize(self, means, std_devs, projection_matrix) -> np.array:
				"""Standadize dataset:
                param: means: np.array
                param: std_devs: np.array
                param: projection_matrix: np.array
                output: reduced_dataset"""
				dataset = self.model_raw_data
				std_dataset = np.zeros(dataset.shape)
				std_dataset = (dataset - means) / std_devs
				reduced_dataset = std_dataset.dot(projection_matrix)
				return reduced_dataset
	def provide_raw_loops(self, transmem_doms_limit):
    
		if transmem_doms_limit:
        
        
			new_seq = ''
			i_lengths = []
			for t, label in zip(self.sequence, self.tmhmm_seq):
				if label in ['o', 'O']:
					new_seq += t
				elif label in ['m', 'M']:
					new_seq += 'X'
				elif label in ['i', 'I']:
					new_seq += ''
				i_lengths.append(1)
                
			avg_i_length = sum(i_lengths) / len(i_lengths) if i_lengths else 0
      
        
			new_seq = new_seq.replace('X', '')
        
            
			return new_seq
		else:
			conds = ['o', 'O']
			if self.localization == "OuterMembrane":
				conds += ['i', 'I']
			new_seq = ""
			for i in range(self.length):
				if self.tmhmm_seq[i] in conds:
					new_seq += self.sequence[i]
				elif len(new_seq) > 0 and not new_seq[len(new_seq)-1] == "X":
					new_seq += "X"
			return new_seq.split('X')

	
	def provide_raw_loops_std(self):
		#print("Warning: this method uses X as a exclusive symbol to split the final protein. Check if X is used inside the protein sequence!")
		conds = ['o', 'O']
		if self.localization == "out":
			conds += ['i', 'I']
		new_seq = ""
		for i in range(self.length):
			if self.tmhmm_seq[i] in conds:
				new_seq += self.sequence[i]
			elif len(new_seq) > 0 and not new_seq[len(new_seq)-1] == "X":
				new_seq += "X"
		return new_seq.split('X')
		
    
	@staticmethod
	def hsp_match_parser(hsp_match, query, parsing_window_size=9, max_sub=3, max_mismatch=1):
		to_return = []
		if max_sub >= parsing_window_size or max_mismatch >= parsing_window_size:
			print("Error in the match parser. Max substitutions and max mismatches have to be lower than the parsing window size! Return.")
			return to_return
		for i in range(0, len(hsp_match)-(parsing_window_size-1)):
			tmp = hsp_match[i:i+parsing_window_size]
			tmp_sub = 0
			tmp_mismatch = 0
			for el in tmp:
				if el == " ":
					tmp_mismatch += 1
				elif el == "+":
					tmp_sub += 1
			if tmp_sub <= max_sub and tmp_mismatch <= max_mismatch:
				to_return.append({'match': tmp, 'query': query, 'start_pos': i})
		return to_return
	
	@staticmethod
	def peptide_comparison(parser_output_seq, peptide):
		match = parser_output_seq['match']
		query = parser_output_seq['query']
		starting_position = parser_output_seq['start_pos']
		tmp_match = query[starting_position:starting_position+len(match)]
		extended_matches = []
		for i in range(len(peptide)-len(tmp_match)+1):
			if tmp_match == peptide[i:i+len(tmp_match)]:
				start_pep, start_query, len_tmp_match = Protein.extendLeft(peptide, i, query, starting_position, len(tmp_match))
				len_tmp_match = Protein.extendRight(peptide, start_pep, query, start_query, len_tmp_match)
				extended_matches.append(query[start_query:start_query+len_tmp_match])
		return extended_matches
	
	@staticmethod
	def extendLeft(peptide, start_pep, real_query, start_q, len_query):
		while (peptide[start_pep:start_pep+len_query] == real_query[start_q:start_q+len_query]) and (start_pep > 0) and (peptide[start_pep-1:start_pep+len_query] == real_query[start_q-1:start_q+len_query]):
			start_pep -= 1
			start_q -= 1
			len_query += 1
		return start_pep, start_q, len_query
	
	@staticmethod
	def extendRight(peptide, start_pep, real_query, start_q, len_query):
		while (peptide[start_pep:start_pep+len_query] == real_query[start_q:start_q+len_query]) and (start_pep+len_query<len(peptide)) and (peptide[start_pep:start_pep+len_query+1] == real_query[start_q:start_q+len_query+1]):
			len_query += 1
		return len_query	
		
	@staticmethod
	def information_to_csv(list_of_proteins):
		print("information_to_csv method call: this method returns a .csv file with the available information about the given proteins, in the same order as in the list.")
		from pandas import DataFrame
		DataFrame([[str(protein.id),
                 str("".join([str(protein.accession) if protein.accession!=None else ""])),
                 str(protein.length),
                 str(protein.transmembrane_doms),
                 str(protein.localization[0].localization),
                 str(protein.localization[0].reliability),
                 #str(", ".join([str(element) for element in protein.localization])),
                 str("".join([str(round(protein.p_vir,4)) if protein.p_vir!=None else ""])),
                 str("".join([str(round(protein.p_ad, 4)) if protein.p_ad!=None else ""])),
                 str("".join([str(round(protein.conservation_score, 4)) if protein.conservation_score!=None else ""])),
                 str("".join(str(len([str(dic['match']) for dic in protein.list_of_shared_human_peps if len(protein.list_of_shared_human_peps)>0])))),
                 str("".join(str(len([str(dic['match']) for dic in protein.list_of_shared_mouse_peps if len(protein.list_of_shared_mouse_peps)>0])))),
                 str("".join(str(len([str(dic['match']) for dic in protein.list_of_shared_conserv_proteome_peps if len(protein.list_of_shared_conserv_proteome_peps)>0])))),
                 str("".join([str(round(protein.sapiens_peptides_sum,4)) if protein.sapiens_peptides_sum!=None else "0"])),
                 str("".join([str(round(protein.mouse_peptides_sum,4)) if protein.mouse_peptides_sum!=None else "0"])),
                 str("".join([str(protein.annotations) if protein.annotations!=None else ""])),
                 str(", ".join(list(set(protein.list_of_peptides_from_comparison_with_mhcpep_sapiens)))), 
                 str(", ".join(list(set(protein.list_of_peptides_from_comparison_with_mhcpep_mouse)))),  
                 str(protein.sequence),
                 str("".join([str(protein.original_sequence_if_razor) if protein.original_sequence_if_razor!=None else ""])),
	         str("".join([str(protein.sequence_out) if protein.sequence_out!=None else ""])),
                 str("".join([str(protein.tmhmm_seq) if "M" in str(protein.tmhmm_seq) else ""])),
		 str("".join([str(protein.MHC1_binders) if protein.MHC1_binders != None else ""])),
                 str("".join([str(protein.MHC2_binders) if protein.MHC2_binders != None else ""])),
                    
                 str("".join([str(protein.MHC1_pb_binders) if protein.MHC1_pb_binders != None else ""])),
                 str("".join([str(protein.MHC2_pb_binders) if protein.MHC2_pb_binders != None else ""])),
                 
                 ] for protein in list_of_proteins
                ], 
                columns= ['id ',
                    'uniprot_accession_code',
                    'length',
                    'transmembrane_doms',
                    'localization',
                    'localization score',
                    'virulence_probability',
                    'adhesin_probability',
                    'conservation_score',
                    'list_of_shared_human_peps',
                    'list_of_shared_mouse_peps',
                    'list_of_shared_conserv_proteome_peps',
                    'human_peptides_sum',
                    'mouse_peptides_sum',
                    'annotations',
                    'list_of_peptides_from_comparison_with_mhcpep_sapiens',
                    'list_of_peptides_from_comparison_with_mhcpep_mouse',
                    'sequence',
                    'original_sequence_if_razor',
		    'sequence_out',
                    'tmhmm_seq'
		    'MHC1_binders',
                    'MHC2_binders',
                    'MHC1_pb_binders',
                    'MHC2_pb_binders',
                   
                     ]
                ).to_csv(outfile) 
