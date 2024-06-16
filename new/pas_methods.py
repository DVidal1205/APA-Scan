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

def Get_Peak_Positions(chromosomes, ann_df, p1_dir, p2_dir, p1_name, p2_name, output_dir, extended):
	print("Getting list of all 3'-end-seq peaks")
	output_columns = ['Chrom', 'Gene', 'Strand', 'Start', 'End',  'Positions']
	writer_list = []
	filename = os.path.join(output_dir, 'Peak_positions.csv')

	for chrom in chromosomes:
		tt = time.time()
		
		bam_list1, position_row1 = methods.read_bamfiles(p1_dir, p1_name, chrom)
		bam_list2, position_row2 = methods.read_bamfiles(p2_dir, p2_name, chrom)

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
						for i in range(exonCount):
							exonStart = int(exonStartList[i])
							exonEnd = int(exonEndList[i])
							if strand == '+':
								if(exonStart<=cdsEnd and cdsEnd<=exonEnd):
									if (exonStart, txEnd) not in targetList:
										targetList.append((exonStart, txEnd))
							elif strand == '-':
								if(exonStart<=cdsStart and cdsStart<=exonEnd):
									if (txStart, exonEnd) not in targetList:
										targetList.append((txStart, exonEnd))

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

				for (st, en) in revisedTargetList:
					(p1, r1) = methods.makeSplittedList(position_row1, bam_list1, st, en)
					(p2, r2) = methods.makeSplittedList(position_row2, bam_list2, st, en)
					listOfPeakRange1 = methods.findPeakPosition(p1, r1)
					listOfPeakRange2 = methods.findPeakPosition(p2, r2)

					if len(listOfPeakRange1) > 0 or len(listOfPeakRange2) > 0:
						cleavageSites = methods.mergePeaksFromBothSamples(listOfPeakRange1, listOfPeakRange2, strand)
						writer_list.append((chrom, gene, strand, st, en, cleavageSites))

		print("Chrom", chrom, "done in", round((time.time() - tt)/60, 2), "minutes")

	df_output = pd.DataFrame(writer_list, columns=output_columns)
	df_output.to_csv(filename, sep='\t')

#done
def with_PA_peaks(chromosomes, s1_dir, s2_dir, g1_name, g2_name, output_dir, result_filename):
	df_p = pd.read_csv(os.path.join(output_dir, "Peak_positions.csv"), delimiter="\t")

	writer_list = []
	output_columns = ['Chrom', 'Gene Name', 'Strand', 'Start', 'End', 'Position', 'p-value', 'Ratio Difference', 'Absolute ratio difference', 'n1: '+g1_name, 'n2: '+g2_name, 'N1: '+g1_name, 'N2: '+g2_name]
	
	position_row = []
	for chrom in chromosomes:
		ss = time.time()
		bam_list1, position_row1 = methods.read_bamfiles(s1_dir, g1_name, chrom)
		bam_list2, position_row2 = methods.read_bamfiles(s2_dir, g2_name, chrom)
		
		target_tt = df_p.loc[df_p['Chrom']==chrom]
		for index, row in target_tt.iterrows():
			gene = row['Gene'].strip()
			strand = row['Strand']
			start = int(row['Start'])
			end = int(row['End'])
			positionsList = list(map(int, row['Positions'].strip('[\' ]').split(',')))
			
			signi_p = 1.1
			ratio_diff = 'Nan'
			signi_ratio_diff = 'Nan'
			abs_ratio_diff = 'Nan'
			n1_f, N1_f, n2_f, N2_f = 0, 0, 0, 0
			targetRC1_f, targetRC2_f, RC1_f, RC2_f = 0, 0, 0, 0
			length_f, targetLength_f, signi_pos = 0, 0, 0
			flag = 0
			for pos in positionsList:
				length = end - start + 1
				targetLength = length
				
				if (strand == '+' and (pos-start)<=(length*0.85)) or (strand == '-' and (pos-start)>=(length*0.15)):
					if strand == '+':
						length = pos-start
						targetLength = end - pos
						RC1 = methods.CountReadCoverage(chrom, start, pos, bam_list1, position_row1)
						RC2 = methods.CountReadCoverage(chrom, start, pos, bam_list2, position_row2)
						targetRC1 = methods.CountReadCoverage(chrom, pos+1, end, bam_list1, position_row1)
						targetRC2 = methods.CountReadCoverage(chrom, pos+1, end, bam_list2, position_row2)
					else:
						targetLength = pos-start
						length = end - pos
						RC1 = methods.CountReadCoverage(chrom, pos+1, end, bam_list1, position_row1)
						RC2 = methods.CountReadCoverage(chrom, pos+1, end, bam_list2, position_row2)
						targetRC1 = methods.CountReadCoverage(chrom, start, pos, bam_list1, position_row1)
						targetRC2 = methods.CountReadCoverage(chrom, start, pos, bam_list2, position_row2)
					
					#print(targetRC1, targetRC2, targetLength, RC1, RC2, length)
					
					n1 = targetRC1/targetLength
					n2 = targetRC2/targetLength
					N1 = RC1/length
					N2 = RC2/length

					#print(n1, n2, N2, N2)

					if N1!=0 and N2!=0:
						ratio_diff = (n1/N1) - (n2/N2)
					else:
						ratio_diff = 0

					N1 = N1 + n1
					N2 = N2 + n2

					if N1!= 0 and N2!=0:
						P0 = (n1+n2)/(N1+N2)
						n10 = N1 * P0
						n20 = N2 * P0
						exp = [n10, N1-n10, n20, N2-n20]
						if 0 not in exp:
							flag = 1
							res = chisquare([n1, N1-n1, n2, N2-n2], f_exp=exp, ddof = 1)
							if res[1] < signi_p:
								signi_p = res[1]
								signi_ratio_diff = ratio_diff
								abs_ratio_diff = abs(signi_ratio_diff)
								n1_f = n1
								N1_f = N1
								n2_f = n2
								N2_f = N2
								RC1_f = RC1
								RC2_f = RC2
								targetRC1_f = targetRC1
								targetRC2_f = targetRC2
								signi_pos = pos
								length_f = length
								targetLength_f = targetLength

			if flag == 1:			
				writer_list.append((chrom, gene, strand, start, end, signi_pos, signi_p, signi_ratio_diff, abs_ratio_diff, n1_f, n2_f, N1_f-n1_f, N2_f-n2_f))
			
		print("Chrom ", chrom, "done in ", round((time.time() - ss)/60, 2), "minutes")
		#sys.exit()

	df_output = pd.DataFrame(writer_list, columns=output_columns)
	df_output.to_csv(os.path.join(output_dir, result_filename+".csv"), sep='\t')

	print("APA-Scan quantification done.")
	return

#done 
def with_PA_peaks_all(chromosomes, s1_dir, s2_dir, g1_name, g2_name, output_dir, result_filename):
	df_p = pd.read_csv(os.path.join(output_dir, "Peak_positions.csv"), delimiter="\t")

	writer_list = []
	output_columns = ['Chrom', 'Gene Name', 'Strand', 'Start', 'End', 'Position', 'p-value', 'Ratio Difference', 'Absolute ratio difference', 'n1: '+g1_name, 'n2: '+g2_name, 'N1: '+g1_name, 'N2: '+g2_name]
	
	position_row = []
	for chrom in chromosomes:
		print(chrom)
		ss = time.time()
		bam_list1, position_row1 = methods.read_bamfiles(s1_dir, g1_name, chrom)
		bam_list2, position_row2 = methods.read_bamfiles(s2_dir, g2_name, chrom)
		
		target_tt = df_p.loc[df_p['Chrom']==chrom]
		for index, row in target_tt.iterrows():
			gene = row['Gene'].strip()
			strand = row['Strand']
			start = int(row['Start'])
			end = int(row['End'])
			positionsList = list(map(int, row['Positions'].strip('[\' ]').split(',')))

			for pos in positionsList:
				pos = int(pos.strip())
				length = end - start + 1
				targetLength = length
				
				if (strand == '+' and (pos-start)<=(length*0.85)) or (strand == '-' and (pos-start)>=(length*0.15)):
					if strand == '+':
						length = pos-start
						targetLength = end - pos
						RC1 = methods.CountReadCoverage(chrom, start, pos, bam_list1, position_row1)
						RC2 = methods.CountReadCoverage(chrom, start, pos, bam_list2, position_row2)
						targetRC1 = methods.CountReadCoverage(chrom, pos+1, end, bam_list1, position_row1)
						targetRC2 = methods.CountReadCoverage(chrom, pos+1, end, bam_list2, position_row2)
					else:
						targetLength = pos-start
						length = end - pos
						RC1 = methods.CountReadCoverage(chrom, pos+1, end, bam_list1, position_row1)
						RC2 = methods.CountReadCoverage(chrom, pos+1, end, bam_list2, position_row2)
						targetRC1 = methods.CountReadCoverage(chrom, start, pos, bam_list1, position_row1)
						targetRC2 = methods.CountReadCoverage(chrom, start, pos, bam_list2, position_row2)
					
					n1 = targetRC1/targetLength
					n2 = targetRC2/targetLength
					N1 = RC1/length
					N2 = RC2/length

					if N1!=0 and N2!=0:
						ratio_diff = (n1/N1) - (n2/N2)
					else:
						ratio_diff = 0

					N1 = N1 + n1
					N2 = N2 + n2

					if N1!= 0 and N2!=0:
						P0 = (n1+n2)/(N1+N2)
						n10 = N1 * P0
						n20 = N2 * P0
						exp = [n10, N1-n10, n20, N2-n20]
						if 0 not in exp:
							flag = 1
							res = chisquare([n1, N1-n1, n2, N2-n2], f_exp=exp, ddof = 1)
							writer_list.append((chrom, gene, strand, start, end, pos, res[1], ratio_diff, abs(ratio_diff), n1, n2, N1-n1, N2-n2))

		print("Chrom", chrom, "done in ", roudn((time.time() - ss)/60, 2), "minutes")
	
	df_output = pd.DataFrame(writer_list, columns=output_columns)
	df_output.to_csv(result_filename, sep='\t')
	print("APA-Scan quantification done.")
	return

def findPeakPosition(p, r):
	listOfPeakPositions = []
	listOfPeakRange = []
	flag = 0
	for i in range(len(p)):
		if int(r[i])>1:
			if flag == 0:
				peakStartPoint = p[i]
				peakAreaPos = []
				peakAreaRead = []
				flag = 1
			if flag == 1:
				peakAreaPos.append(p[i])
				peakAreaRead.append(r[i])

		if (int(r[i])==1 and flag == 1) or (i == len(p)-1 and flag == 1):
			if i == len(p)-1:
				peakEndPoint = p[i]
			else:
				peakEndPoint = p[i-1]

			np_peakAreaPos = np.array(peakAreaPos)
			np_peakAreaRead = np.array(peakAreaRead)

			indexes = peakutils.indexes(np_peakAreaRead, thres=6, min_dist=35, thres_abs = True)

			numberOfIndex = len(indexes)

			if numberOfIndex == 1:
				listOfPeakRange.append((peakAreaPos[indexes[0]], peakStartPoint, peakEndPoint))

			elif numberOfIndex > 1:
				peakRangeList = calculatePeaksWithVallys(peakAreaPos, peakAreaRead, indexes, peakStartPoint, peakEndPoint)
			
				listOfPeakRange.extend(peakRangeList)

			flag = 0

	return listOfPeakRange

def mergePeaksFromBothSamples(listOfPeakRange1, listOfPeakRange2, strand):
	len1 = len(listOfPeakRange1)
	len2 = len(listOfPeakRange2)

	cleavageSites = []
	peakFlag = []
	minVal = 10000

	if len1 == 0:
		for (peak, st, en) in listOfPeakRange2:
			if strand == '+':
				cleavageSites.append(en)
			else:
				cleavageSites.append(st)
		return cleavageSites
	elif len2 == 0:
		for (peak, st, en) in listOfPeakRange1:
			if strand == '+':
				cleavageSites.append(en)
			else:
				cleavageSites.append(st)
		return cleavageSites

	if len1<len2:
		row = len1
		col = len2
		firstList = listOfPeakRange1
		secondList = listOfPeakRange2
	else:
		row = len2
		col = len1
		firstList = listOfPeakRange2
		secondList = listOfPeakRange1

	peakFlag1 = {}
	peakFlag2 = {}
	distanceToSortList = []

	for i in range(row):
		for j in range(col):
			if strand == '+':
				distance = abs(int(firstList[i][2])-int(secondList[j][2]))
				distanceToSortList.append((int(firstList[i][2]), int(secondList[j][2]), distance))
				peakFlag1[firstList[i][2]] = 0
				peakFlag2[secondList[j][2]] = 0
			else:
				distance = abs(int(firstList[i][1])-int(secondList[j][1]))
				distanceToSortList.append((int(firstList[i][1]), int(secondList[j][1]), distance))
				peakFlag1[firstList[i][1]] = 0
				peakFlag2[secondList[j][1]] = 0

	distanceToSortList.sort(key = itemgetter(2), reverse = False)
	
	for (pos1, pos2, dist) in distanceToSortList:
		if peakFlag1[pos1] == 0 and peakFlag2[pos2] == 0:
			averageOfTwoPeaks = math.ceil((pos1+pos2)/2)
			cleavageSites.append(averageOfTwoPeaks)
			peakFlag1[pos1] = 1
			peakFlag2[pos2] = 1
		elif peakFlag1[pos1] != 0 and peakFlag2[pos2] == 0:
			cleavageSites.append(pos2)
			peakFlag2[pos2] = 1
		elif peakFlag1[pos1] == 0 and peakFlag2[pos2] != 0:
			cleavageSites.append(pos1)
			peakFlag1[pos1] = 1

	return cleavageSites
