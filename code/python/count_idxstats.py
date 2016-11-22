#!/usr/bin/env python 

# USAGE: count_idxstats.py input_file read_length
# Counts number of mapped reads in an idxstats file and computes the 90% cutoff

import sys

with open(sys.argv[1], 'r') as infile:
	
	total = 0 

	for line in infile:
	
		if line.startswith('*'): continue
	
		line = line.split()
		total += int(line[1])

pick = int(round(float(total) * 0.9))
print('\nTotal reads: ' + str(total))	
print('90% of total: ' + str(pick) + '\n')
