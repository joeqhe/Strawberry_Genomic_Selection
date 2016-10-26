#!/usr/bin/python
#SNP_cluster_identification.py
#input required: csv files with SNP locations as integer number
#identifies number of SNPs in an n length fragment at positions 1 -> n
#repeats for 2 -> (n+1) ... until end of sequence
#outputs top x sections

n = 1000
x = 200
file = "vesca_snp_positions.csv"

import csv

##########################################################################################
#Reads data file as dictionary
##########################################################################################

row_count = 0
LG_count = {}
SNP_data = {}

with open(file) as csvfile:
	readCSV = csv.reader(csvfile, delimiter = ",")
	
	for row in readCSV:
		SNP_id, LG, position = row
		SNP_data [SNP_id] = [int(LG),int(position)]
		if SNP_data [SNP_id] [0] not in LG_count:
			LG_count [SNP_data [SNP_id] [0]] = 1
		else:
			LG_count [SNP_data [SNP_id] [0]] = LG_count [SNP_data [SNP_id] [0]] + 1
		
		row_count = row_count + 1
		
print "this file contains %d SNPs" % (row_count)

for row in LG_count:
	print "with %d SNPs in LG %d" % (LG_count [int(row)], int(row))

##########################################################################################
#
##########################################################################################

SNP_locations = []

for SNP_name in SNP_data:
	SNP_locations.append (SNP_data [SNP_name] [1])
max_SNP_pos = max(SNP_locations)

print "the highest position of any SNP is %d" % (max_SNP_pos)

##########################################################################################
#Counts number of SNPs in windows of size n
##########################################################################################


SNP_count = {}

for window_start in range (1,max_SNP_pos):
	print "working on %d of %d windows" % (window_start, max_SNP_pos)
	for SNP_name in SNP_data:
		for LG in range (1,len(LG_count)):
			if window_start <= SNP_data [SNP_name][1] <= (window_start + n):
				if "LG %d, start %d, end %d" % (LG, window_start, window_start +n) not in SNP_count:
					SNP_count ["LG %d, start %d, end %d" % (LG, window_start, window_start +n)] = 1
				else:
					SNP_count ["LG %d, start %d, end %d" % (LG, window_start, window_start +n)] = SNP_count ["LG %d, start %d, end %d" % (LG, window_start, window_start +n)] + 1





