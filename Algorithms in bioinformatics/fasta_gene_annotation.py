#!/usr/bin/env python3

"""
Author: Thierry Haddad
Date: 23-11-2018
Description: Annotates potential genes for a given FASTA file.
			 For every gene, GC percentage is calculated and stop codon
			 is presented. Additionally, the average GC percentage per
			 stop codon is calculated.
"""

from sys import argv, exit
from subprocess import check_call
import os.path

def parse_fasta():
	""" Returns dictionary of FASTA file in format."""
    # Get file argument, catch wrong/missing entries
	try:
		input_file = argv[1]
		# Make sure file is FASTA format
		assert input_file.endswith('.fasta') or input_file.endswith('.fa')
	except IndexError:
		print('No file argument given')
		exit(1)
	except AssertionError:
		print('Given file is not in FASTA format')
		exit(1)
	except Exception as e:
		print(str(e))
		exit(1)

	# Get FASTA sequences
	with open(input_file, 'r') as inp:
		lines = inp.readlines()

	# Dictionary for the FASTA content
	fasta = {}
	id_count = 0
	for line in lines:
		# Don't add empty lines
		if not line.strip():
			continue
		# Header line
		if line.startswith('>'):
			id_count += 1
			fasta[id_count] = {}
			header = line.strip()[1:]
			fasta[id_count]['header'] = header
			fasta[id_count]['sequence'] = ""
		# Sequence line
		else:
			fasta[id_count]['sequence'] += line.strip()
	return fasta, input_file

def run_augustus(input_file=''):
	"""
	Function to run Augustus on the given FASTA file.
	Augustus will perform annotation on the FASTA sequences,
	giving predictions for potential protein-coding genes.

	input_file: String name of the FASTA file
	"""

	# Only run if output file doesn't already exist
	if not os.path.isfile('aug.gff'):
		command = ("augustus --strand=both --genemodel=complete --gff3=on "
				"--species=aspergillus_nidulans "
				"{}".format(input_file))
		# Run command and save as GFF3 file
		with open('aug.gff','w') as ao:
			check_call(command, stdout=ao, shell=True)
	return True

def parse_gff():
	"""
	Opens the .gff file created by Augustus and parses the content to a
	Dictionary. Saves the gene number, sequence header, start location,
	stop location, intron locations and strand.
	"""

	# Input file
	if os.path.isfile('aug.gff'):
		with open('aug.gff', 'r') as aug:
			lines = aug.readlines()
	else:
		print('GFF file missing! Exiting.')
		exit(1)
	
	# Dict for parsed .gff content
	gff = {}

	# Go through all GFF lines
	for line in lines:
		# Empty line
		if not line.strip():
			continue
		# Catch header for FASTA sequences
		elif 'name = ' in line:
			header = line.rstrip().split('name = ')[-1].split(')')[0]
		# Start of a new gene
		elif line.startswith('# start gene g'):
			gen_id = line.rstrip().split(' ')[-1]  # Gene number as unique id
			gff[gen_id] = {}
			gff[gen_id]['header'] = header
			gff[gen_id]['intron'] = []
		# Catch start and stop codons
		elif 'start_codon' in line or 'stop_codon' in line:
			line = line.rstrip().split('\t')
			# Start or stop codon
			codon_type = line[2].split('_')[0]
			gff[gen_id][codon_type] = {}
			# Codon locations
			loc1 = int(line[3])
			loc2 = int(line[4])
			gff[gen_id][codon_type]['locations'] = [loc1, loc2]
			# Strand (+/-) of codon
			gff[gen_id][codon_type]['strand'] = line[6]
		# Catch introns
		elif 'intron' in line:
			gff[gen_id]['intron'].append(line.split('\t')[3:5])
	return gff

def get_sequences(gff={}, fasta={}):
	"""
	Retrieve nucleotide sequences from the FASTA file based on the locations
	and strand found per gene in the gff Dictionary.

	gff: Dictionary from Augustus containing:
	per gene id:
	 - location 1
	 - location 2
	 - strand

	fasta: Dictionary from FASTA file containing:
	per sequence:
	 - header
	 - nucleotide sequence
	"""

	# Nucleotide dictionary for reverse complement
	nuc_dict = {
			'A':'T', 'T':'A',
			'C':'G', 'G':'C',
			}
	# Go through FASTA sequences
	for seq_name in fasta.keys():
		# Keep GFF file in order
		for i in range(1, len(gff.keys())+1):
			gene = 'g{}'.format(i)
			# Direction of gene
			strand = gff[gene]['start']['strand']
			start = gff[gene]['start']['locations'][0 if strand=='+' else 1]
			stop = gff[gene]['stop']['locations'][1 if strand=='+' else 0]
			header = gff[gene]['header']
			# Only search if gene originates from current FASTA seq
			if header in fasta[seq_name]['header']:
				# Forward strand
				if strand == '+':
					seq = fasta[seq_name]['sequence'][start-1:stop]
					stop_codon = seq[-3:]
				# Reverse strand
				else:
					seq = fasta[seq_name]['sequence'][stop-1:start]
					pre_stop_codon = seq[0:3]
					stop_codon = ''
					for nuc in pre_stop_codon:
						nuc = nuc_dict[nuc]
						stop_codon += nuc
					stop_codon = stop_codon[::-1]
				# Remove introns
				for intron in gff[gene]['intron']:
					loc1 = int(intron[0])-1
					loc2 = int(intron[1])
					int_seq = seq[loc1:loc2]
					# Cut specific intron region
					seq = "".join(seq.split('int_seq'))
				# Calculate GC percentage
				gc = float(sum(seq.count(nuc) for nuc in ['C','G']))
				gc_percent = gc/len(seq)*100
				gff[gene]['gc'] = gc_percent
				gff[gene]['stopcodon'] = stop_codon
	return gff

def write_tab_del(gff={}):
	"""
	Prints a tab-delimited output for genes in the GFF Dict.
	For all genes it prints the gene number, gc percentage and stop codon.
	
	gff: Dictionary containing the annotated genes.
	"""

	# Dict for stop codon gc percentage
	stops = {}
	for i in range(1, len(gff.keys())+1):
		gene = 'g{}'.format(i)
		gc = '{0:.3f}'.format(gff[gene]['gc'])
		stop = gff[gene]['stopcodon']
		if stop not in stops:
			stops[stop] = [float(gc)]
		else:
			stops[stop].append(float(gc))
		line = '{}\t{}\t{}'.format(gene, gc, stop)
		print(line)
	print('---')
	
	# Average GC percentage per stop codon type
	for stop in stops.keys():
		avg_gc = '{0:.3f}'.format(sum(stops[stop])/len(stops[stop]))
		print('{}\t{}'.format(stop, avg_gc))
	return None

if __name__ == '__main__':
	# Get FASTA content
	fasta, input_file = parse_fasta()
	# Run Augustus for gene annotation
	if run_augustus(input_file):
		# Parse GFF contents from Augustus
		gff = parse_gff()
		# Get sequence characteristics
		gff = get_sequences(gff, fasta)
		# Print output in tab-delimited format
		write_tab_del(gff)
