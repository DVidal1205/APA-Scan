import csv
import time
import sys
from operator import itemgetter
import math
import pandas as pd
from Bio import SeqIO
import re
import bisect
from bisect import bisect_left
from scipy.stats import chisquare
import numpy as np
import methods
import rna_methods
import pas_methods
import glob, os
import configparser

# Set flags
(speciesFlag, inputFlag, outFlag, pasFlag) = (0, 0, 0, 0)

# Set the start time
startTime = time.time()

# Function to list directories
def list_dirs(path):
    return [os.path.basename(x) for x in filter(os.path.isdir, glob.glob(os.path.join(path, '*')))]

# Read the configuration file
config = configparser.ConfigParser()
config.read('configuration.ini')

# Parse the configuration file and handle the values
input1_dir = config['INPUT_RNAseq']['input1']
input2_dir = config['INPUT_RNAseq']['input2']
if input1_dir[-1] == "/":
	input1_dir = input1_dir[:-1]
if input2_dir[-1] == "/":
	input2_dir = input2_dir[:-1]
pasSeq_dir1 = config['INPUT_PASseq']['pas1']
pasSeq_dir2 = config['INPUT_PASseq']['pas2']
if pasSeq_dir1[:-1] == "/":
	pasSeq_dir1 = pasSeq_dir1[:-1]
if pasSeq_dir2[:-1] == "/":
	pasSeq_dir2 = pasSeq_dir2[:-1]
output_dir = config['OUTPUT_FOLDER']['output_dir']
if output_dir[-1] != "/":
	output_dir += "/"
extended = config['Extended_3UTR']['extended']
all_events = config['All_events']['All']

# Create the output directory
os.makedirs(output_dir, exist_ok=True)

# Set the annotation, reference genome and input directories
inp_annotation = config['ANNOTATION']['annotation']
ref_genome = config['ANNOTATION']['genome']
g1_name, g2_name = os.path.basename(input1_dir), os.path.basename(input2_dir)

# Print the configuration values for the user
print("RNA-seq input 1 dir:", input1_dir)
print("RNA-seq input 2 dir:", input2_dir)
print("3'-end-seq input 1 dir:", pasSeq_dir1)
print("3'-end-seq input 2 dir:", pasSeq_dir2)
print("Output Dir:", output_dir) 
print("Annotation:", inp_annotation, ref_genome, "\n")
print("Loading annotation information ...")

# Load the annotation file and extract the chromosomes
ann_df = pd.read_csv(inp_annotation, delimiter='\t')
chr_list = ann_df['chrom'].str.split("_", n = 1, expand = True)
chromosomes = list(set(chr_list[0])-set(['chrM', 'chrUn']))

# Clean the annotation file
if '#name' in ann_df.columns:
	ann_df.rename(columns = {"#name": "name"}, inplace = True)

# Generate the read coverage files
readTime = time.time()
print("Generating read coverage files from the RNA-seq data...")
for sample1 in os.listdir(input1_dir):
	if sample1.endswith('.bam'):
		methods.SamtoText(input1_dir, sample1, chromosomes)
for sample2 in os.listdir(input2_dir):
	if sample2.endswith('.bam'):
		methods.SamtoText(input2_dir, sample2, chromosomes)
print("Read coverage files generated in", round((time.time() - readTime)/60, 2), "minutes.")

########################################################################
# BEGIN APA-SCAN PIPELINE										       #
########################################################################

# Create the result filename
result_filename = "APA_Scan_"+g1_name+"_Vs_"+g2_name

# Check if the 3'-end-seq data is available
if pasSeq_dir1 == 'NULL' or pasSeq_dir2=='NULL':
	s1_namelist = list_dirs(input1_dir)
	s2_namelist = list_dirs(input2_dir)
	
	print("Preparing result using RNA-seq data only")
	rnaTime = time.time()
	rna_methods.Get_Signal_Positions(chromosomes, ann_df, ref_genome, output_dir, extended)
	if all_events == 'all':
		rna_methods.with_PAS_signal_all(chromosomes, input1_dir, input2_dir, s1_namelist, s2_namelist, g1_name, g2_name, output_dir, result_filename)
	else:
		rna_methods.with_PAS_signal(chromosomes, input1_dir, input2_dir, s1_namelist, s2_namelist, g1_name, g2_name, output_dir, result_filename)
	print("Total time (without PAS seq data):", round((time.time() - rnaTime)/60, 2), "minutes.")

else:
	print("Generating read coverage files for the 3'-end-seq data...")
	rnaPasTime = time.time()
	for sample1 in os.listdir(pasSeq_dir1):
		if sample1.endswith('.bam'):
			methods.SamtoText(pasSeq_dir1, sample1, chromosomes)

	for sample2 in os.listdir(pasSeq_dir2):
		if sample2.endswith('.bam'):
			methods.SamtoText(pasSeq_dir2, sample2, chromosomes)
	
	p1_namelist = list_dirs(pasSeq_dir1)
	p2_namelist = list_dirs(pasSeq_dir2)
	p1_name = pasSeq_dir1.split("/")[-1]
	p2_name = pasSeq_dir2.split("/")[-1]

	pas_methods.Get_Peak_Positions(chromosomes, ann_df, pasSeq_dir1, pasSeq_dir2, p1_name, p2_name, output_dir, extended)
	if all_events == 'all':
		pas_methods.with_PA_peaks_all(chromosomes, input1_dir, input2_dir, g1_name, g2_name, output_dir, result_filename)
	else:
		pas_methods.with_PA_peaks(chromosomes, input1_dir, input2_dir, g1_name, g2_name, output_dir, result_filename)
	print("Total time (with PAS seq data):", round((time.time() - rnaPasTime)/60, 2), "minutes.")

print("Total time:", round((time.time() - startTime)/60, 2), "minutes.")
read_file = pd.read_csv(os.path.join(output_dir, result_filename+".csv"), delimiter = '\t')
read_file.to_excel (os.path.join(output_dir, result_filename+".xlsx"), index = None, header=True)

