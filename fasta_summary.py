#!/usr/bin/python

# generates a summary of fasta files
# lists first 10 characters in sequence along with length and identification
# change line 8 to target appropriate file

from Bio import SeqIO
handle = open("/home/hejoe/GT_seq_primer_design/a.lines.fasta","rU")
for record in SeqIO.parse(handle,"fasta"):
	print(record.id)
	print("amplicon length:",len(record))
	print(record[:10])
handle.close()


