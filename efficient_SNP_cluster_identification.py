#!/usr/bin/python
#input required: csv files with SNP locations at linkage group and integer position
#identifies number of SNPs in an n length fragment for all fragments containing >0 SNPs
#keeps first fragment with highest number of SNPs and deletes those SNPs from SNP list
#repeats for x times for densest SNP locations
#outputs top x sections

n = 100000
x = 100
file = "vesca_snp_positions.csv"

import csv
import copy

result = {}

##########################################################################################
#Full scan function
##########################################################################################


def full_scan (SNP_data, bins, n):
	best_window = [None, None, 0]
	#best_window = [window_start, SNP_id, SNP_count]
	for bin_id in bins:
		SNP_set = bins [bin_id]
		if bin_id + 1 in bins:
			SNP_set = SNP_set | bins [bin_id +1]
		for SNP_id in bins [bin_id]:
			window_start = SNP_data [SNP_id]
			current_window = [window_start,[],0]
			for SNP in SNP_set:
				if window_start <= SNP_data [SNP] <= window_start + n:
					current_window [1].append (SNP)
					current_window [2] += 1
				if best_window [2] < current_window [2]:
					best_window = current_window

	result [best_window[0]] = [best_window [1],best_window [2]]
	
	for SNP_id in best_window [1]:
		while SNP_id in bins [SNP_data [SNP_id] / n]:
			bins [SNP_data [SNP_id] / n].remove (SNP_id)
		
		del SNP_data [SNP_id]

	
##########################################################################################
#Reads data file as dictionary
##########################################################################################

row_count = 0
LG_count = {}
SNP_data = {}

with open(file) as csvfile:
	readCSV = csv.reader(csvfile, delimiter = ",")
	
	for row_count,row in enumerate(readCSV):
		SNP_id, LG, position = row
		LG = int (LG)
		SNP_data [SNP_id] = LG * 100000000 + int(position)
#		SNP_data key in form (LG)12345678 where the numbers is the SNP position
		if LG not in LG_count:
			LG_count [LG] = 1
		else:
			LG_count [LG] +=  1
		
SNP_data_original = copy.deepcopy (SNP_data)

print "this file contains %d SNPs" % (row_count)

for LG in LG_count:
	print "with %d SNPs in LG %d" % (LG_count [LG], LG)
	
	
##########################################################################################
#Create bins
##########################################################################################

bins = {}

for SNP_id in SNP_data:
	bin_id  = SNP_data [SNP_id] / n
	if bin_id not in bins:
		bins [bin_id] = set([SNP_id])
	else:
		bins [bin_id].add(SNP_id)

		
##########################################################################################
#Scan entire genome
##########################################################################################


for results in range (1,x+1):
	print "calculating %d of %d densest %d nucleotides in %s" % (results, x, n, file)
	full_scan (SNP_data, bins, n)
	
print "this file contains %d SNPs" % (row_count)

for LG in LG_count:
	print "with %d SNPs in LG %d" % (LG_count [LG], LG)


##########################################################################################
#Remove overlapping SNP regions
##########################################################################################

for entry_count, entry in enumerate (result):
	print "filtering %d of %d results for overlap" % (entry_count, len(result))
	overlap_regions = {}
	remove_counter = 0
	for key in result.keys():
		if (entry - n) <= key <= (entry + n):
			overlap_regions [key] = result [key]
						
	if len (overlap_regions) > 1:
		density = []
		
		for overlap in overlap_regions:
			density.append (overlap_regions [overlap][-1])
			
		densest = max (density)
		
		for overlap in overlap_regions:
			if overlap_regions [overlap][-1] == densest:
				keep_SNP = overlap
				
		for key in overlap_regions:
			if key is not keep_SNP:
				del result [key]
				remove_counter += 1
				
print "%d result(s) removed due to overlap" % (remove_counter)
print "type 'result' to see densest SNP regions"	



