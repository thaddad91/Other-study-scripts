#!/usr/bin/env python3

"""
Author: Thierry Haddad
Student nr: 911031296110
Description: Parses a given FASTQ file, calculates sequences characteristics,
			 Phred scores and per-position read qualities. Also trims the 
			 FASTQ file and compares the per-position scores to the original
			 qualities.
"""

from sys import argv, exit
import os.path
from subprocess import check_output

def get_input(input_file=''):
	"""
	Opens a given .fq file, parses content to a Dictionary and returns it.
	The following is saved: 
	 - header, starting with '@'
	 - sequence
	 - '+', optional characteristics
	 - per-position quality score

	input file: String that contains the file name.
	"""
	try:
		# Go through the FASTQ file line by line
		with open(input_file, "r") as inp:
			fq = {}
			i = 0
			header = ''  # Starts with '@'
			for line in inp:
				if i == 0 or not i % 4:  # 4 Lines per sequence
					if line.startswith('@'):
						header = line.rstrip('\n')
						fq[header] = []
				else:
					if header != '':  # Skip empty lines
						fq[header].append(line.rstrip('\n'))
				i += 1
		return fq

	except IndexError:
		print('Please add the .fq file as argument.')
		exit(1)

	except Exception as e:
		print(str(e))
		exit(1)

def calc_phred_quality(fq={}):
	"""
	Calculates per-position Phred score of the quality line, added as
	a list of Integers in the fq Dictionary.

	fq: FASTQ Dictionary containing the following:
	 - header starting with '@'
	 - sequence
	 - '+' followed by optional characteristics
	 - per-position quality scores
	"""
	for header in fq:
		quality_line = list(fq[header][-1])
		phred = []
		for score in quality_line:
			# Illumina 1.5 = ASCII+64, so detracting 64
			p_score = ord(score)-64  
			if p_score > 41:
				print('Invalid Phred score found (symbol/score):')
				print(score, p_score)
				exit(1)
			phred.append(p_score)
		fq[header].append(phred)
	return fq

def trim_fq(input_file = ''):
	"""
	Runs command line tool 'fastq_quality_trimmer' to trim the 
	given FASTQ file.

	input_file: String that contains the name of the given FASTQ file.
	"""
	try:
		if os.path.isfile(input_file):  # Check if file exists
			command = [
				'fastq_quality_trimmer',
				'-t', '30',  # Threshold 30
				'-Q', '64',  # ASCII+64
				'-i', input_file,  # Input file
				'-o', 'trimmed.fq',  # Output file
			]		
			res = check_output(command)
			return True

	except Exception as e:
		print(str(e))
		exit(1)

def calc_seq_length(fq={}):
	"""
	Calculates the minimum, maximum and average read length
	for both FASTQ files. Appends it to the FASTQ Dictionary.

	fq: Dictionary of the FASTQ file containing the following:
	 - header
	 - sequence
	 - optional characteristics
	 - per-position quality scores
	"""
	seq_lengths = []
	# FASTQ header starting with '@'
	for header in fq:
		# Sequence length
		seq_lengths.append(len(fq[header][0]))
	max_seq = max(seq_lengths)
	min_seq = min(seq_lengths)
	avg_seq = sum(seq_lengths)/float(len(seq_lengths))
	fq["lengths"] = [min_seq, max_seq, avg_seq]
	return fq

def calc_score_per_pos(fq={}, trimmed_fq={}, longest_read=0):
	"""
	Compared the per-position scores of both FASTQ files.
	Notes the original average, trimmed average and net difference.

	fq: Dictionary of the FASTQ file:
	 - header
	 - sequence
	 - optinal characteristics
	 - per-position quality scores
	trimmed_fq: Identical to Dict fq, but for the trimmed file.
	longest_read: Integer of length of the longest read.
	"""
	positions = longest_read  # Maximum amount of positions
	for pos in range(positions):  # pos = index position in read
		# Original FASTQ
		ori_score = []
		for i in fq.keys():  # For every header
			try:
				ori_score.append(fq[i][3][pos])  # Append position pos
			except IndexError:  # Shorter reads
				pass
		avg_ori = sum(ori_score)/float(len(ori_score))

		# Trimmed FASTQ
		trim_score = []
		for i in trimmed_fq.keys():
			try:
				trim_score.append(trimmed_fq[i][3][pos])
			except IndexError:  # Shorter reads
				pass
		avg_trim = sum(trim_score)/float(len(trim_score))

		diff = avg_trim-avg_ori
		nrs = [pos+1, avg_ori, avg_trim, diff]
		print('{}\t{}\t{}\t{}'.format(*nrs))
	return None

def print_seq_lengths(tf={}, trimmed_fq={}):
	"""
	"""
	fq_len = fq["lengths"]
	trimmed_fq_len = trimmed_fq["lengths"]
	print('ORIGINAL: min={}, max={}, avg={}'.format(*fq_len))
	print('TRIMMED: min={}, max={}, avg={}'.format(*trimmed_fq_len))
	longest_read = fq["lengths"][1]
	del fq["lengths"]
	del trimmed_fq["lengths"]
	return fq, trimmed_fq, longest_read

if __name__ == '__main__':
	try:
		input_file = argv[1]
		assert input_file.endswith('.fq')
	except IndexError:  # Catch no arguments
		print('Please add FASTQ file as argument 1')
		exit(1)
	except AssertionError:  # Catch wrong file type
		print('Argument file is not a .fq file')
		exit(1)
	except Exception as e:  # Other errors
		print(str(e))
		exit(1)
	# Normal FASTQ file
	fq = get_input(input_file)
	fq = calc_phred_quality(fq)
	# Trimmed FASTQ file
	trim_fq(input_file)
	trimmed_fq = get_input('trimmed.fq')
	trimmed_fq = calc_phred_quality(trimmed_fq)
	# Compare read lengths
	fq = calc_seq_length(fq)
	trimmed_fq = calc_seq_length(trimmed_fq)
	fq, trimmed_fq, longest_read = print_seq_lengths(fq, trimmed_fq)
	# Compare per-position scores
	calc_score_per_pos(fq, trimmed_fq, longest_read)


