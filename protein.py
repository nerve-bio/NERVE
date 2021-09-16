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
		self.conservation = None
		self.function = None
		self.scan_prosite_information = None
		
	def print_information(self):
		print("Information about protein " + str(self.id) + ":")
		print("   accession number = " + str(self.accession))
		print("   length = " + str(self.length))
		print("   localization = " + str(self.localization))
		print("   estimated probability to be an adhesin = " + str(self.p_ad))
		print("   number of transmembrane domains = " + str(self.transmembrane_doms))
		print("   number of peptides shared with sapiens = " + str(self.sapiens_peptides_sum))
		print("   number of peptides shared with mhcpep = " + str(self.conservation) +"\n")
		print("   putative function of the protein = not yet implemented\n")
	
	@staticmethod
	def hsp_match_parser(hsp_match, parsing_window_size=9, max_sub=3, max_mismatch=1):
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
				to_return.append(tmp)
		return to_return
