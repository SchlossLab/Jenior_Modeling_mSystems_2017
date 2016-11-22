#!/usr/bin/env python 

# USAGE: normalize_idxstats.py input_file read_length
# Normalizes mapped read count to library and target sizes
# Ignores discordant read mappings

import sys

read = int(sys.argv[2])

with open(sys.argv[1], 'r') as infile:
	
	out_str = str(sys.argv[1]).rstrip('txt') + 'norm.txt'
	outfile = open(out_str, 'w')
	
	for line in infile:
	
		if line.startswith('*'): continue
	
		line = line.split()
		target = str(line[0])
		length = int(line[1])
		mapped = int(line[2])
		
		normalized = (mapped * read) / length
		
		entry = target + '\t' + str(normalized) + '\n'
		outfile.write(entry)
		
outfile.close()
