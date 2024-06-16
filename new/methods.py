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

class Stack:
	def __init__(self):
		self.items = []

	def size(self):
		return len(self.items)

	def isEmpty(self):
		return self.items == []

	def push(self, val):
		self.items.append(val)

	def top(self):
		if self.isEmpty():
			return None
		else:
			return self.items[self.size()-1]

	def pop(self):
		if self.isEmpty():
			return None
		else:
			return self.items.pop()

def bi_contains(lst, item):
    return bisect_left(lst, item)

def SamtoText(input_path, bamfile, chromosomes):
	output_dir = os.path.join(input_path, bamfile[:-4])
	samtools_dir = "/usr/bin/samtools-0.1.8/samtools"
	os.makedirs(output_dir, exist_ok=True)
	cwd = os.getcwd()
	cmd1 = samtools_dir+" index "+os.path.join(input_path, bamfile)		# make samtools index bamfile.bam.bai
	os.system(cmd1)

	print(bamfile,"...")
	for chrom in chromosomes:
		cmd2 = samtools_dir+" view -b "+os.path.join(input_path, bamfile)+" "+chrom+" -o "+os.path.join(output_dir, chrom+".bam")
		cmd3 = samtools_dir+" pileup "+os.path.join(output_dir, chrom+".bam")+" | cut -f 2,4 > "+os.path.join(output_dir, chrom+".txt")   ### Need to use pileup, not mpileup
		command = cmd2+";"+cmd3
		os.system(command)
	return

def getFinalTargetRegion(inputList):
	n = len(inputList)
	inputList.sort(key = itemgetter(1), reverse = True)
	st = Stack()
	st.push(inputList[0])

	for i in range(1,n):
		stacktop = st.top()
		if inputList[i][0] < stacktop[0] and inputList[i][1] < stacktop[0]:
			st.push(inputList[i])
		elif inputList[i][0] < stacktop[0] and inputList[i][1] > stacktop[0] and inputList[i][1] < stacktop[1]:
			st.pop()
			st.push((inputList[i][0], stacktop[1]))
		elif inputList[i][1] == stacktop[1] and inputList[i][0] < stacktop[0]:
			st.pop()
			st.push(inputList[i])

	newList = []
	while(True):
		if st.size() == 0:
			break;
		stacktop = st.top()
		newList.append(stacktop)
		st.pop()

	return newList


def findAllOccurance(p, s):
    i = s.find(p)
    while i != -1:
        yield i
        i = s.find(p, i+1)

def CountReadCoverage(chrom, start, end, bam_list, position_row):
	totalCount = 0

	pos1 = bi_contains(position_row, start)
	pos2 = bi_contains(position_row, end)

	if pos2 >= len(bam_list) or int(bam_list[pos2][0]) != end:
		pos2 = pos2 - 1

	if(pos1 < len(bam_list) and pos2 < len(bam_list)):
		for t in range(pos1, pos2+1):
			read = int(bam_list[t][1])
			totalCount += read
		
	return totalCount

def read_bamfiles(input_dir, sample, chrom):
	ss = time.time()
	bam_df = pd.read_csv(os.path.join(input_dir, sample, chrom+".txt"), delimiter='\t')
	position_row = bam_df.iloc[:, 0].tolist()
	bam_list = bam_df.values.tolist()
	return bam_list, position_row

def Generate_coverage(chrom, start, end, pos, strand, bam_list, position_row):
	if strand == '+':
		length = pos-start
		targetLength = end-pos
		RC = CountReadCoverage(chrom, start, pos, bam_list, position_row)
		targetRC = CountReadCoverage(chrom, pos+1, end, bam_list, position_row)
	else:
		targetLength = pos-start
		length = end-pos
		RC = CountReadCoverage(chrom, pos+1, end, bam_list, position_row)
		targetRC = CountReadCoverage(chrom, start, pos, bam_list, position_row)
	
	n = targetRC/targetLength
	N = RC/length
	return n, N

def makeSplittedList(position_row, bam_list, start, end):						
	p = []
	r = []

	pos1 = bi_contains(position_row, start)
	pos2 = bi_contains(position_row, end)

	if pos2 >= len(bam_list) or int(bam_list[pos2][1]) != end:
		pos2 = pos2 - 1

	if(pos1 < len(bam_list) and pos2 < len(bam_list)):
		for t in range(pos1, pos2+1):
			p.append(int(bam_list[t][0]))
			r.append(int(bam_list[t][1]))

	return (p,r)

def calculatePeaksWithVallys(peakAreaPos, peakAreaRead, indexes, allStartPoint, allEndPoint):
	peakRangeList = []
	peakStartPoint = allStartPoint
	prevPeakValue = 0
	indexLen = len(indexes)

	# because we will compare two peaks together to consider the valley in between, we wont take the last peak index in the for loop
	for k in range(indexLen-1):
		selectedRange = peakAreaRead[indexes[k]:indexes[k+1]+1]
		valley = min(selectedRange)
		valleyPos = selectedRange.index(valley)
		actualValleyPos = peakAreaPos[indexes[k]+valleyPos]

		if valley < 0.3*min(peakAreaRead[indexes[k]],peakAreaRead[indexes[k+1]]):
			peakEndPoint = actualValleyPos
			peakRangeList.append((peakAreaPos[peakAreaRead.index(max(prevPeakValue,peakAreaRead[indexes[k]]))], peakStartPoint, peakEndPoint))
			prevPeakValue = 0
			peakStartPoint = actualValleyPos + 1
		else:
			prevPeakValue = max(prevPeakValue, peakAreaRead[indexes[k]])


	ind = [i for i, value in enumerate(peakAreaRead) if value == max(prevPeakValue, peakAreaRead[indexes[indexLen-1]]) and peakAreaPos[i] >= peakStartPoint]

	# because ind is a list with indexes of the maxpeak value, for that we will take ind[0] (ind[0] or ind[1] are the indices of the same value)
	peakRangeList.append((peakAreaPos[ind[0]], peakStartPoint, allEndPoint))
	return peakRangeList


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


