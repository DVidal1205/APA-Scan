import csv
import time
import sys
import math
from operator import itemgetter
import pandas as pd
from Bio import SeqIO
import re
import bisect
from bisect import bisect_left
from scipy.stats import chisquare
import peakutils
import numpy as np
import os
import methods

def Get_Signal_Positions(chromosomes, ann_df, ref_genome, output_dir, extended):
    filename = os.path.join(output_dir, 'Signal_positions.csv')
    writer_list = []
    output_columns = ['Chrom', 'Gene', 'Strand', 'Start', 'End', 'Positions']
    fasta_sequences = SeqIO.parse(open(ref_genome),'fasta')

    for fasta in fasta_sequences:
    	chrom, sequence = fasta.id, str(fasta.seq)
    	if chrom in chromosomes:
    		tt = time.time()
    		geneList = list(set(ann_df[ann_df['chrom']==chrom]['name2'].tolist()))

    		for gene in geneList:
    			gene_rows = ann_df.loc[(ann_df['chrom'] == chrom) & (ann_df['name2'] == gene)]
    			targetList = []
    			for index, row in gene_rows.iterrows():
    				if str(row['name'].strip()).startswith('NM_') or str(row['name'].strip()).startswith('NR_'):
    					cdsStart, cdsEnd = int(row['cdsStart']), int(row['cdsEnd'])
    					txStart, txEnd = int(row['txStart']), int(row['txEnd'])
    					exonCount, strand = int(row['exonCount']), row['strand']
    					exonStartList = row['exonStarts'].split(',')[0:-1]
    					exonEndList = row['exonEnds'].split(',')[0:-1]

    					if cdsStart < cdsEnd:
	    					if strand == '+':
	    						regStart = cdsEnd
	    						regEnd = txEnd
	    						for i in range(exonCount):
	    							exonStart = int(exonStartList[i])
	    							exonEnd = int(exonEndList[i])
	    							if(exonStart<=regStart and regStart<=exonEnd):
	    								if (exonStart, regEnd) not in targetList:
	    									targetList.append((exonStart, regEnd))
	    									break

	    					elif strand == '-':
	    						regStart = txStart
	    						regEnd = cdsStart
	    						for i in range(exonCount):
	    							exonStart = int(exonStartList[i])
	    							exonEnd = int(exonEndList[i])
	    							if(exonStart<=regEnd and regEnd<=exonEnd):
	    								if (regStart, exonEnd) not in targetList:
	    									targetList.append((regStart, exonEnd))
	    									break

	    		if len(targetList) > 0:				
    				revisedTargetList = methods.getFinalTargetRegion(targetList)

    				if extended.upper() == 'YES':
    					if strand == '+':
    						st, en = revisedTargetList[-1]
    						del revisedTargetList[-1]
    						revisedTargetList.append((st, en+10000))
    					elif strand == '-':
    						st, en = revisedTargetList[0]
    						del revisedTargetList[0]
    						revisedTargetList.append((st-10000, en))

    				if strand == '+':
    					for (st, en) in revisedTargetList:
    						length = en - st + 1
    						seq = sequence[st-2: en-2].upper()
    						listPlus = []
    						for i in methods.findAllOccurance('AATAAA', seq):
    							pos = i+5
    							if pos<=(length*0.85):
    								listPlus.append(st+pos)
    						for i in methods.findAllOccurance('ATTAAA', seq):
    							pos = i+5
    							if pos<=(length*0.85):
    								listPlus.append(st+pos)

    						if len(listPlus)>0:
    							writer_list.append((chrom, gene, strand, st, en, list(set(listPlus))))

    				elif strand == '-':
    					for (st, en) in revisedTargetList:
    						length = en - st + 1
    						seq = sequence[st-2: en-2].upper()
    						listMinus = []
    						for pos in methods.findAllOccurance('TTTATT', seq):
    							if pos>=(length*0.15):
    								listMinus.append(st+pos)
    						for pos in methods.findAllOccurance('TTTAAT', seq):
    							if pos>=(length*0.15):
    								listMinus.append(st+pos)

    						if len(listMinus)>0:
    							writer_list.append((chrom, gene, strand, st, en, list(set(listMinus))))

    		print("Chrom",chrom, " done in ", round((time.time() - tt)/60, 2), "minutes")

    df_output = pd.DataFrame(writer_list, columns=output_columns)
    df_output.to_csv(filename, sep="\t")
    
def with_PAS_signal(chromosomes, input1_dir, input2_dir, s1_namelist, s2_namelist, g1, g2, output_dir, result_filename):
	len1,len2 = len(s1_namelist), len(s2_namelist)

	df_p = pd.read_csv(os.path.join(output_dir, "Signal_positions.csv"), delimiter="\t")
	writer_list = []
	output_columns = ['Chrom', 'Gene', 'Strand', 'Start', 'End', 'Position', 'p-value', 'Ratio Difference', 'Absolute ratio difference', 'n1: '+g1, 'n2: '+g2, 'N1: '+g1, 'N2: '+g2]

	position_row = []
	for chrom in chromosomes:
		ss = time.time()
		s1_bam_list, s2_bam_list, s1_position_row, s2_position_row = {}, {}, {}, {}
		for sample1 in s1_namelist:
			bam_list, position_row = methods.read_bamfiles(input1_dir, sample1, chrom)
			s1_bam_list[sample1] = bam_list
			s1_position_row[sample1] = position_row
		for sample2 in s2_namelist:
			bam_list, position_row = methods.read_bamfiles(input2_dir, sample2, chrom)
			s2_bam_list[sample2] = bam_list
			s2_position_row[sample2] = position_row

		selected_rows = df_p.loc[df_p['Chrom']==chrom]
		for index, row in selected_rows.iterrows():
			gene = row['Gene'].strip()
			strand = row['Strand']
			start = int(row['Start'])
			end = int(row['End'])
			positionsList = list(map(int, row['Positions'].strip('[\' ]').split(',')))
			#print("positionsList", positionsList)
			signi_p = 1.1
			ratio_diff = 'Nan'
			signi_ratio_diff = 'Nan'
			abs_ratio_diff = 'Nan'
			n1_f, N1_f, n2_f, N2_f, signi_pos = 0, 0, 0, 0, 0
			flag = 0

			for pos in positionsList:
				length = end - start + 1
				targetLength = length
				
				s1_n, s1_N, s2_n, s2_N = 0, 0, 0, 0
				for sample1 in s1_namelist:
					n, N = methods.Generate_coverage(chrom, start, end, pos, strand, s1_bam_list[sample1], s1_position_row[sample1])
					s1_n += n
					s1_N += N
				for sample2 in s2_namelist:
					n, N = methods.Generate_coverage(chrom, start, end, pos, strand, s2_bam_list[sample2], s2_position_row[sample2])
					s2_n += n
					s2_N += N

				n1 = s1_n/len1
				N1 = s1_N/len1
				n2 = s2_n/len2
				N2 = s2_N/len2
				#print("n1, N1, n2, N2", n1, N1, n2, N2)
				
				ratio_diff = 0

				N1 = N1 + n1
				N2 = N2 + n2

				if N1>0 and N2>0:
					ratio_diff = (n1/N1) - (n2/N2)
					P0 = (n1+n2)/(N1+N2)
					n10 = N1 * P0
					n20 = N2 * P0
					exp = [n10, N1-n10, n20, N2-n20]
					if 0 not in exp:
						flag = 1
						res, p_value = chisquare([n1, N1-n1, n2, N2-n2], f_exp=exp, ddof = 1)
						#print("res, p_value", res, p_value)
						if p_value < signi_p:
							signi_p = p_value
							signi_ratio_diff = ratio_diff
							abs_ratio_diff = abs(signi_ratio_diff)
							n1_f = n1
							N1_f = N1
							n2_f = n2
							N2_f = N2
							signi_pos = pos
				
			if flag == 1:
				writer_list.append((chrom, gene, strand, start, end, signi_pos, signi_p, signi_ratio_diff, abs_ratio_diff, n1_f, n2_f, N1_f-n1_f, N2_f-n2_f))
		print("Chrom ", chrom, " done in ", round((time.time() - ss)/60, 2), "minutes")

	df_output = pd.DataFrame(writer_list, columns=output_columns)
	df_output.to_csv(os.path.join(output_dir, result_filename+".csv"), sep='\t')
	print("APA-Scan quantification done.")
	return

def with_PAS_signal_all(chromosomes, input1_dir, input2_dir, s1_namelist, s2_namelist, g1, g2, output_dir, result_filename):
	len1, len2 = len(s1_namelist), len(s2_namelist)

	df_p = pd.read_csv(os.path.join(output_dir, "Signal_positions.csv"), delimiter="\t")
	writer_list = []
	output_columns = ['Chrom', 'Gene', 'Strand', 'Start', 'End', 'Position', 'p-value', 'Ratio Difference', 'Absolute ratio difference', 'n1: '+g1, 'n2: '+g2, 'N1: '+g1, 'N2: '+g2]

	position_row = []
	for chrom in chromosomes:
		ss = time.time()
		s1_bam_list, s2_bam_list, s1_position_row, s2_position_row = {}, {}, {}, {}
		for sample1 in s1_namelist:
			bam_list, position_row = methods.read_bamfiles(input1_dir, sample1, chrom)
			s1_bam_list[sample1] = bam_list
			s1_position_row[sample1] = position_row
		for sample2 in s2_namelist:
			bam_list, position_row = methods.read_bamfiles(input2_dir, sample2, chrom)
			s2_bam_list[sample2] = bam_list
			s2_position_row[sample2] = position_row

		selected_rows = df_p.loc[df_p['Chrom']==chrom]
		for index, row in selected_rows.iterrows():
			gene = row['Gene'].strip()
			strand = row['Strand']
			start = int(row['Start'])
			end = int(row['End'])
			positionsList = list(map(int, row['Positions'].strip('[\' ]').split(',')))

			for pos in positionsList:
				pos = int(pos.strip())
				length = end - start + 1
				targetLength = length
				
				s1_n, s1_N, s2_n, s2_N = 0, 0, 0, 0
				for sample1 in s1_namelist:
					n, N = methods.Generate_coverage(chrom, start, end, pos, strand, s1_bam_list[sample1], s1_position_row[sample1])
					s1_n += n
					s1_N += N
				for sample2 in s2_namelist:
					n, N = methods.Generate_coverage(chrom, start, end, pos, strand, s2_bam_list[sample2], s2_position_row[sample2])
					s2_n += n
					s2_N += N

				n1 = s1_n/len1
				N1 = s1_N/len1
				n2 = s2_n/len2
				N2 = s2_N/len2
				
				ratio_diff = 0

				N1 = N1 + n1
				N2 = N2 + n2

				if N1>0 and N2>0:
					ratio_diff = (n1/N1) - (n2/N2)
					P0 = (n1+n2)/(N1+N2)
					n10 = N1 * P0
					n20 = N2 * P0
					exp = [n10, N1-n10, n20, N2-n20]
					if 0 not in exp:
						flag = 1
						res, p_value = chisquare([n1, N1-n1, n2, N2-n2], f_exp=exp, ddof = 1)
						writer_list.append((chrom, gene, strand, start, end, pos, p_value, ratio_diff, abs(ratio_diff), n1, n2, N1-n1, N2-n2))

		print("Chrom", chrom, "done in ",  round((time.time() - ss)/60, 2), "minutes")

	df_output = pd.DataFrame(writer_list, columns=output_columns)
	df_output.to_csv(os.path.join(output_dir, result_filename+".csv"), sep='\t')
	print("APA-Scan quantification done.")
	return
