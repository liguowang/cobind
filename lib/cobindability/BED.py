#!/usr/bin/env python

import sys
#from bx.bitset import *
from bx.bitset_builders import binned_bitsets_from_file, binned_bitsets_from_list
#from bx.intervals import *
from cobindability import ireader
import logging
import numpy as np
import pandas as pd

__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="0.1.0"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Development"



def unionBed3(lst):
	'''Take the union of 3 column bed files. return a new list'''
	bitsets = binned_bitsets_from_list(lst)
	ret_lst=[]
	for chrom in bitsets:
		bits = bitsets[chrom]
		end = 0
		while 1:
			start = bits.next_set( end )
			if start == bits.size: break
			end = bits.next_clear( start )
			ret_lst.append([chrom, start, end])
	bitsets=dict()
	return ret_lst

def intersectBed3(lst1,lst2):
	'''Take the intersection of two bed files (3 column bed files)'''
	bits1 = binned_bitsets_from_list(lst1)
	bits2 = binned_bitsets_from_list(lst2)

	bitsets = dict()
	ret_lst = []
	for key in bits1:
		if key in bits2:
			bits1[key].iand( bits2[key] )
			bitsets[key] = bits1[key]

	for chrom in bitsets:
		bits = bitsets[chrom]
		end = 0
		while 1:
			start = bits.next_set( end )
			if start == bits.size: break
			end = bits.next_clear( start )
			ret_lst.append([chrom, start, end])
	bits1.clear()
	bits2.clear()
	bitsets.clear()
	return ret_lst

def subtractBed3(lst1,lst2):
	'''subtrack lst2 from lst1'''
	bitsets1 = binned_bitsets_from_list(lst1)
	bitsets2 = binned_bitsets_from_list(lst2)

	ret_lst=[]
	for chrom in bitsets1:
		if chrom not in bitsets1:
			continue
		bits1 = bitsets1[chrom]
		if chrom in bitsets2:
			bits2 = bitsets2[chrom]
			bits2.invert()
			bits1.iand( bits2 )
		end=0
		while 1:
			start = bits1.next_set( end )
			if start == bits1.size: break
			end = bits1.next_clear( start )
			ret_lst.append([chrom,start,end])
	bitsets1 = dict()
	bitsets2 = dict()
	return ret_lst

def tillingBed(chrName,chrSize,stepSize=10000):
	'''tilling whome genome into small sizes'''
	#tilling genome
	for start in range(0,chrSize,stepSize):
		end = start + stepSize
		if end < chrSize:
			yield (chrName,start,end)
		else:
			yield (chrName,start,chrSize)

def bed_actual_size(bed1, bed2):
	'''
	Calculate the total size of BED file.

	Examples
	--------
	$cat test.bed
		chr1  0  100
		chr1  50  150
		chr1  80  180

	Note:
	1. The *actual* size of 'test.bed' is 300 (bp). However, the *genomic* size is
	180 (bp), because the overlapped possitions were only counted once!
	2. "bedfile1" and "bedfile2" can be regular, compressed, or
	remote files.
	3. "bedfile1" and "bedfile2" can be bigBed files. The suffix for bigBed file
	must be one of ('.bb','.bigbed','.bigBed','.BigBed', '.BB',' BIGBED').

	Parameters
	----------
	bed1 : str
		File name of the first BED file.
	bed2 : str
		File name of the second BED file.

	Returns
	-------
	List (bed1_size, bed2_size, bed1_count, bed2_count).

	'''
	bed_file1_size = 0
	bed_file2_size = 0
	bed_file1_count = 0
	bed_file2_count = 0
	#print >>sys.stderr, "reading %s ..." %	 bedfile1
	for l in ireader.reader(bed1):
		l = l.strip()
		if l.startswith(('browser','#','track')):
			continue
		f = l.split()
		if (int(f[2]) - int(f[1])) < 0:
			logging.error("invalid BED line: %s" % l)
			sys.exit(1)
		bed_file1_size += (int(f[2]) - int(f[1]))
		bed_file1_count += 1
	#print >>sys.stderr, "reading %s ..." %	 bedfile2
	for l in ireader.reader(bed2):
		l = l.strip()
		if l.startswith(('browser','#','track')):
			continue
		f = l.split()
		if (int(f[2]) - int(f[1])) < 0:
			logging.error("invalid BED line: %s" % l)
			sys.exit(1)
		bed_file2_size += (int(f[2]) - int(f[1]))
		bed_file2_count += 1
	return	(bed_file1_size, bed_file2_size, bed_file1_count, bed_file2_count)

def bed_genomic_size(*argv):
	'''
	Calculate the *genomic* size of BED file.

	Examples
	--------
	$cat test.bed
		chr1  0  100
		chr1  50  150
		chr1  80  180
	Note:
	1. The *actual* size of 'test.bed' is 300 (bp). However, the *genomic* size is
	180 (bp), because the overlapped possitions were only counted once!
	2. "bedfile1" and "bedfile2" can be regular, compressed, or
	remote files.
	3. "bedfile1" and "bedfile2" can be bigBed files. The suffix for bigBed file
	must be one of ('.bb','.bigbed','.bigBed','.BigBed', '.BB',' BIGBED').
	Parameters
	----------
	argv : list of genomic regions.
		Each argument can be bed, bed-like or bigBed format file.


	Returns
	-------
	Pandas Series.
	'''
	union_sizes = []
	for arg in argv:
		if type(arg) is list:
			bitsets = binned_bitsets_from_list(arg )
		else:
			bitsets = binned_bitsets_from_file( ireader.reader(arg) )
		union_size = 0
		for chrom in bitsets:
			bits = bitsets[chrom]
			end = 0
			while 1:
				start = bits.next_set( end )
				if start == bits.size: break
				end = bits.next_clear( start )
				union_size += (end - start)
		union_sizes.append(union_size)
	return (union_sizes)

def bedinfo(infile):
	logging.debug("Gathering teh basic statistics of BED file: %s" % infile)
	bed_infor={}
	bed_infor['Name'] = infile
	bed_infor['Genomic_size'] = bed_genomic_size(infile)[0]
	bed_infor['Total_size'] = 0
	bed_infor['Count'] = 0
	size = 0
	sizes = []
	for l in ireader.reader(infile):
		if l.startswith(('browser','#','track')):
			continue
		f = l.split()
		if len(f) < 3:
			logging.error("invalid BED line: %s" % l)
		if (int(f[2]) - int(f[1])) < 0:
			logging.error("invalid BED line: %s" % l)
		bed_infor['Count'] += 1
		size = int(f[2]) - int(f[1])
		bed_infor['Total_size'] += size
		sizes.append(size)
	bed_infor['Mean_size'] = np.mean(sizes)
	bed_infor['Median_size'] = np.median(sizes)
	bed_infor['Min_size'] = np.min(sizes)
	bed_infor['Max_size'] = np.max(sizes)

	return(pd.Series(data = bed_infor, index=['Name','Count', 'Genomic_size', 'Total_size', 'Max_size', 'Mean_size', 'Median_size', 'Min_size']))

def bed_overlap_size(bed1,bed2):
	"""
	Calculate the total number of *bases* overlapped between two bed files

	Parameters
	----------
	bed1 : str or list
		File name of the first BED file. Can also be a list, such as
		[(chr1 100 200), (chr2 150  300), (chr2 1000 1200)]
	bed2 : str or list
		File name of the second BED file. Can also be a list, such as
		[(chr1 100 200), (chr2 150  300), (chr2 1000 1200)]

	Returns
	-------
	Int. Overlapped size.

	"""
	overlap_size = 0
	if type(bed1) is list:
		bits1 = binned_bitsets_from_list( bed1 )
	else:
		bits1 = binned_bitsets_from_file( ireader.reader(bed1) )
	if type(bed2) is list:
		bits2 = binned_bitsets_from_list( bed2 )
	else:
		bits2 = binned_bitsets_from_file( ireader.reader(bed2) )

	bitsets = dict()
	if bed1 == bed2:
		return 0.0
	for key in bits1:
		if key in bits2:
			bits1[key].iand( bits2[key] )
			bitsets[key] = bits1[key]

	for chrom in bitsets:
		bits = bitsets[chrom]
		end = 0
		while 1:
			start = bits.next_set( end )
			if start == bits.size: break
			end = bits.next_clear( start )
			overlap_size += end - start
	return overlap_size





def bedtolist(bedfile):
	'''
	Calculate the total size of BED file.
	'''
	regions = []
	#print >>sys.stderr, "reading %s ..." %	 bedfile1
	for l in ireader.reader(bedfile):
		l = l.strip()
		if l.startswith('browser') or l.startswith('#') or l.startswith('track'):
			continue
		f = l.split()
		if len(f) < 3:
			print("invalid BED file:%s" % l, file=sys.stderr)
			sys.exit(1)
		if (int(f[2]) - int(f[1])) < 0:
			print("invalid BED file:%s" % l, file=sys.stderr)
			sys.exit(1)
		regions.append(  (f[0], int(f[1]), int(f[2]))  )
	return regions

