from Bio import SeqIO
from Bio.Alphabet import IUPAC
chp_list = list(SeqIO.parse("srf_chip.fasta", "fasta", IUPAC.unambiguous_dna))

count = 0
for dna in chp_list:
	match = dna.seq.upper().count('GCCCATATATGG') # .upper vs. .upper()
	count = count+match

import re

chp_match = 0
for dna in chp_list:
	if re.search(r'[GT][CA]CC[AT]TATA[AT]GG', str(dna.seq)): # dna vs. str(dna.seq))
		chp_match = chp_match+1

from Bio import motifs

srf_m = motifs.read(open("MA0083.1.sites"), "sites")
srf_m.pseudocounts = 1
srf_m.background = 0.4

for dna in chp_list: # put into for loop
	for pos, score in srf_m.pssm.search(dna.seq, threshold=7.0): # use .seq atrribute
		print "Position %d: score = %5.2f" % (pos, score)
