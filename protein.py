class protein:
	def __init__(self, identifier, sequence_string):
		self.id = identifier
		self.accession = identifier.split('|')[1]
		self.sequence = sequence_string
		self.length = len(sequence_string)
		self.localization = None
		self.p_ad = 0
		self.transmembrane_doms = None
		self.tmhmm_seq = None
		self.sapiens_peptides_sum = None
		self.conservation = None
		
	def print_information(self):
		print("Information about protein " + str(self.id) + ":")
		print("   accession number = " + str(self.accession))
		print("   length = " + str(self.length))
		print("   localization = " + str(self.localization))
		print("   estimated probability to be an adhesin = " + str(self.p_ad))
		print("   number of transmembrane domains = " + str(self.transmembrane_doms))
		print("   number of peptides shared with sapiens = " + str(self.sapiens_peptides_sum))
		print("   number of peptides shared with mhcpep = " + str(self.conservation))