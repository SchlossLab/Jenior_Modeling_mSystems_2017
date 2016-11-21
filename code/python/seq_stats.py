#!/usr/bin/python2.7
'''USAGE: fasta_stats.py seqFile 
This script calculates various statistics about the provided fasta or fastq file.
'''
import sys
import os

# This function reads in fasta file, appends the length of each sequence to a list, and counts all Gs & Cs.
	# It returns a sorted list of sequence lengths with the G+C total as the last element.
def readFasta(FastaFile):
	
	tempseq = ''
	GC_count = 0
	N_count = 0
	lenList = []

	for line in FastaFile: 
	
		if line.startswith('>'): 
			lenList.append(len(tempseq))
			GC_count = GC_count + tempseq.count('C') + tempseq.count('c') + tempseq.count('G') + tempseq.count('g')
			N_count = N_count + tempseq.count('N') + tempseq.count('n')
			tempseq = ''
			continue

		else:
			tempseq = tempseq + line.strip()
	
	lenList.append(len(tempseq))
	
	del lenList[0]
	lenList.sort()
	lenList.append(GC_count)
	lenList.append(N_count)
	
	FastaFile.close()
	
	return lenList


# This function reads in fastq files and returns all sequence lengths in a list
def readFastq(FastqFile):
	
	GC_count = 0
	N_count = 0
	lenList = []
	line_count = 3 

	for line in FastqFile: 
		
		line_count += 1
		
		if line_count == 5:
			seq = line.strip()
			lenList.append(len(seq))
			GC_count = GC_count + seq.count('C') + seq.count('c') + seq.count('G') + seq.count('g')
			N_count = N_count + seq.count('N') + seq.count('n')
			line_count = 1
			continue
		
		else:
			continue
		
	lenList.sort()
	lenList.append(GC_count)
	lenList.append(N_count)
	
	FastqFile.close()
	
	return lenList
	
		
# This function calculates and returns all the printed statistics.
def calcStats(lengths):
	
	Ns = lengths[-1]
	del lengths[-1]
	GCs = lengths[-1] # extract saved GC count
	del lengths[-1] 
	                  	
	total_seq = len(lengths) # Total number of sequences
	len_sum = sum(lengths) # Total number of residues
	total_Mb = len_sum/1000000.00 # Total number of residues expressed in Megabases
	if len_sum == 0: len_sum = 1
	GC_content = (float(GCs)/float(len_sum))*100 # GC content as a percent of total residues

	# interquartile range
	if total_seq >= 4:
		median_len = lengths[int(round(total_seq/2))] # Median sequence length
		q1_range = lengths[0:(int(total_seq/2)-1)]
		if len(q1_range) == 0:
			q1 = lengths[0]
		else:
			q1 = q1_range[int(len(q1_range)/2)]
		q3_range = lengths[(int(total_seq/2)+1):-1]
		if len(q3_range) == 0:
			q3 = lengths[-1]
		else:
			q3 = q3_range[int(len(q3_range)/2)]
		iqr = int(q3 - q1)
	else:
		iqr = 'Too few sequences to calculate'
		median_len = 'Too few sequences to calculate'

	#n50 calculation loop
	current_bases = 0
	n50 = 0
	n90 = 0
	seqs_1000 = 0
	seqs_5000 = 0
	percent50_bases = int(round(len_sum*0.5))
	percent90_bases = int(round(len_sum*0.1))
	
	for object in lengths:
	
		current_bases += object
		
		if object > 1000:
			seqs_1000 += 1
		if object > 5000:
			seqs_5000 += 1
			
		if current_bases >= percent50_bases and n50 == 0:
			n50 = object
		if current_bases >= percent90_bases and n90 == 0:
			n90 = object

	if total_seq < 4:
		n50 = 'Too few sequences to calculate'
		n90 = 'Too few sequences to calculate'
		l50 = 'Too few sequences to calculate'
	else:	
		l50 = lengths.count(n50)

	return(total_seq, total_Mb, n50, median_len, iqr, GC_content, n90, seqs_1000, seqs_5000, Ns, l50)

#----------------------------------------------------------------------------------------#

if os.stat(str(sys.argv[1])).st_size == 0:
	print('ERROR: Empty input file.')
	sys.exit()

input_file = open(sys.argv[1], 'r')

first_line = input_file.readline()

if first_line[0] == '>':
	file_type = 'Fasta'	
	seq_lengths = readFasta(input_file)  
	
elif first_line[0] == '@':	
	file_type = 'Fastq'
	seq_lengths = readFastq(input_file)

else:
	print('ERROR: Invalid file format provided.')
	sys.exit()

stat_list = calcStats(seq_lengths)


output_string = """# Input file name: {filename}
# File type: {filetype}
# Total sequences: {total_seq}
# Total bases: {total_mb} Mb
# Sequence N50: {n50}
# Sequence L50: {l50}
# Sequence N90: {n90}
# Median sequence length: {median_len}
# Interquartile range: {iqr}
# Shortest sequence length: {short}
# Longest sequence length: {long}
# Sequences > 1 kb: {seqs_1000}
# Sequences > 5 kb: {seqs_5000}
# G-C content: {gc}%
# Ns included: {ns}
""".format(filename = str(sys.argv[1]).split('/')[-1], filetype = file_type, total_seq = stat_list[0], total_mb = "%.2f" % stat_list[1], n50 = stat_list[2], median_len = stat_list[3], iqr = stat_list[4], short = seq_lengths[0], long = seq_lengths[-1], gc = "%.2f" % stat_list[5], n90 = stat_list[6], seqs_1000 = stat_list[7], seqs_5000 = stat_list[8], ns = stat_list[9], l50 = stat_list[10])

print output_string
