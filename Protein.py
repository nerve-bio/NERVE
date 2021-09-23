class Protein:
	def __init__(self, identifier, sequence_string):
		self.id = identifier
		self.accession = identifier.split('|')[1]
		self.sequence = sequence_string
		self.length = len(sequence_string)
		self.localization = None
		self.p_ad = 0
		self.transmembrane_doms = None
		self.tmhmm_seq = None
		self.transmembrane_flag = 0 # 1 if transmembrane doms > 2
		self.sapiens_peptides_sum = None
		self.conservation_score = None
		self.function = None
		self.scan_prosite_information = None
		self.list_of_shared_human_peps = []
		
	def print_information(self):
		print("Information about protein " + str(self.id) + ":")
		print("   accession number = " + str(self.accession))
		print("   length = " + str(self.length))
		print("   localization = " + str(self.localization))
		print("   estimated probability to be an adhesin = " + str(self.p_ad))
		print("   number of transmembrane domains = " + str(self.transmembrane_doms))
		if len(self.list_of_shared_human_peps) == 0:
			print("   no interesting peptides shared with sapiens")
		else:
			print("   list of interesting peptides shared with sapiens = " + str(self.list_of_shared_human_peps))
		#print("   number of interesting peptides shared with sapiens = " + str(self.sapiens_peptides_sum))
		print("   conservation with the given proteome blastp score = " + str(self.conservation_score))
		print("   putative function of the protein = not yet implemented\n")
		
	
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
		query = parser_outputs_seq['query']
		starting_position = parser_output_seq['start_pos']
		tmp_match = query[starting_position:starting_position+len(match)]
		extended_matches = []
		for i in range(len(peptide)-len(tmp_match)+1):
			if tmp_match == peptide[i:i+len(tmp_match)]:
				start_pep, start_query, len_tmp_match = extendLeft(peptide, i, query, starting_position, len(tmp_match))
				len_tmp_match = extendRight(peptide, start_pep, query, start_query, len_tmp_match)
				extended_matches.append(query[start_query:start_query+len_tmp_match])
		return list(dict.fromkeys(extended_matches))
	
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
