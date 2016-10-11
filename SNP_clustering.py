#!/usr/bin/python
#input required: csv files with SNP locations as integer number
#identifies number of SNPs in an n length fragment at positions 1 -> n
#repeats for 2 -> (n+1) ... until end of sequence
#outputs top x sections

n = 1000
x = 200


import csv

with open("vesca_snp_positions.csv") as csvfile:
	readCSV = csv.reader(csvfile, delimiter=",")
	for row in readCSV:
		print (row)
	

