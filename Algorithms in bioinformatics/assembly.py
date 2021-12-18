#!/usr/bin/env python3

"""
Author: Thierry Haddad
"""

from sys import argv, exit
from subprocess import check_output, check_call
import os

def assemble_reads(fq_1, fq_2):
	"""x"""
	# Hash paired-end reads
	command1 = [
		"velveth", "Assem", "29", 
		"-shortPaired", "-fastq", "-separate",
		fq_1, fq_2,
	]
	# Make assembly
	command2 = [
		"velvetg", "Assem",
		"-min_contig_lgth", "200",
		"-ins_length", "300",
	]
	try:
		output1 = check_output(command1)
		output2 = check_output(command2)
	except Exception as e:
		print(str(e))
		exit(1)
	return output1, output2

def run_augustus():
	"""x"""
	command = ("augustus --species=saccharomyces_cerevisiae_S288C "
			   "Assem/contigs.fa")
	with open('aut_out','w') as ao:
		check_call(command, stdout=ao, shell=True)

def run_hisat2():
	"""x"""
	command = "hisat2-build yeast/chr3.fa yeast/index/chr3"
	check_call(command, shell=True)
	print("built")
	command = ("hisat2 -p 16 -x yeast/index/chr3 -U "
			"yeast/CENPK_RNA_1.fastq,yeast/CENPK_RNA_2.fastq "
			"-S hisat2_rna.sam")
	check_call(command, shell=True)
	print("mapped rna reads")
	command = ("hisat2 -p 16 -x yeast/index/chr3 -U "
			"yeast/imw004-chr3_1.fastq,yeast/imw004-chr3_2.fastq "
			"-S hisat2_geno.sam")
	check_call(command, shell=True)
	print("mapped genomic reads")

def run_samtools():
	"""FIRST NEEDS .BAM FROM BOWTIE2"""
	# Sorted BAM
	command = "samtools sort -o hisat2_rna.bam hisat2_rna.sam"
	check_call(command, shell=True)
	# Create index file
	command = "samtools index hisat2_rna.bam"
	check_call(command, shell=True)

if __name__ == '__main__':
	try:
		fq_1 = argv[1]
		fq_2 = argv[2]
	except IndexError:  # Catch no arguments
		print('Please add FASTQ file as argument 1')
		exit(1)
	except Exception as e:  # Other errors
		print(str(e))
		exit(1)
	output1, output2 = assemble_reads(fq_1, fq_2)
	#run_augustus()
	#run_hisat2()
	run_samtools()
	print("Done.")
