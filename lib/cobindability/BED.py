#!/usr/bin/env python

import sys,os
import re

import logging
import numpy as np
import pandas as pd
import math
from scipy.stats import fisher_exact
from bx.bitset_builders import binned_bitsets_from_file, binned_bitsets_from_list
from bx.intervals.intersection import Interval, Intersecter
from cobindability import ireader


__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="0.1.0"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Development"



def unionBed3(inbed):
	"""
	Union or merge BED regions.

	Parameters
	----------
	inbed : str or list
		Name of a BED file or list of BED regions, for example,
		[(chr1 100 200), (chr2 150  300), (chr2 1000 1200)]

	Returns
	-------
	ret_lst : list
		List of bed regions with the overlapped regions merged.

	Examples
	--------
	>>> unionBed3([('chr1', 1, 10), ('chr1', 3, 15), ('chr1', 20, 35), ('chr1', 20, 50)])
	[('chr1', 1, 15), ('chr1', 20, 50)]

	"""
	if type(inbed) is list:
		bitsets = binned_bitsets_from_list(inbed)
	elif type(inbed) is str:
		try:
			bitsets = binned_bitsets_from_file( ireader.reader(inbed) )
		except:
			logging.error("invalid input: %s" % inbed)
			sys.exit(1)
	else:
		logging.error("invalid input: %s" % inbed)
		sys.exit(1)
	ret_lst=[]
	for chrom in bitsets:
		bits = bitsets[chrom]
		end = 0
		while 1:
			start = bits.next_set( end )
			if start == bits.size: break
			end = bits.next_clear( start )
			ret_lst.append((chrom, start, end))
	bitsets=dict()
	return ret_lst


def intersectBed3(inbed1,inbed2):
	"""
	Interset two BED files (or lists).

	Parameters
	----------
	inbed1 : str or list
		Name of a BED file or list of BED regions, for example,
		[(chr1 100 200), (chr2 1000 1200)]
	inbed2 : str or list
		Name of a BED file or list of BED regions, for example,
		[(chr1 150 220), (chr2 1100 1300)]

	Returns
	-------
	ret_lst : list
		List of bed regions shared between the two input BED files ( or lists).

	Examples
	--------
	>>> intersectBed3([('chr1', 1, 10), ('chr1', 20, 35)], [('chr1',3, 15), ('chr1',20, 50)])
	[('chr1', 3, 10), ('chr1', 20, 35)]
	"""
	if type(inbed1) is list:
		bits1 = binned_bitsets_from_list(inbed1)
	elif type(inbed1) is str:
		try:
			bits1 = binned_bitsets_from_file(ireader.reader(inbed1))
		except:
			logging.error("invalid input: %s" % inbed1)
			sys.exit(1)
	else:
		logging.error("invalid input: %s" % inbed1)
		sys.exit(1)
	if  type(inbed2) is list:
		bits2 = binned_bitsets_from_list(inbed2)
	elif  type(inbed2) is str:
		try:
			bits2 = binned_bitsets_from_file(ireader.reader(inbed2))
		except:
			logging.error("invalid input: %s" % inbed2)
			sys.exit(1)
	else:
		logging.error("invalid input: %s" % inbed2)
		sys.exit(1)
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
			ret_lst.append((chrom, start, end))
	bits1.clear()
	bits2.clear()
	bitsets.clear()
	return ret_lst


def subtractBed3(inbed1,inbed2):
	"""
	Subtract inbed2 from inbed1.

	Parameters
	----------
	inbed1 : str or list
		Name of a BED file or list of BED regions, for example,
		[(chr1 100 200), (chr2 1000 1200)]
	inbed2 : str or list
		Name of a BED file or list of BED regions, for example,
		[(chr1 150 220), (chr2 1100 1300)]

	Returns
	-------
	ret_lst : list
		List of bed regions of inbed1 with those shared regions with inbed2 removed.

	Examples
	--------
	>>> subtractBed3([('chr1', 1, 10), ('chr1', 20, 35)], [('chr1',3, 15), ('chr1',20, 50)])
	[('chr1', 1, 3)]

	"""
	if type(inbed1) is list:
		bitsets1 = binned_bitsets_from_list(inbed1)
	elif type(inbed1) is str:
		try:
			bitsets1 = binned_bitsets_from_file(ireader.reader(inbed1))
		except:
			logging.error("invalid input: %s" % inbed1)
			sys.exit(1)
	else:
		logging.error("invalid input: %s" % inbed1)
		sys.exit(1)
	if type(inbed2) is list:
		bitsets2 = binned_bitsets_from_list(inbed2)
	elif  type(inbed2) is str:
		try:
			bitsets2 = binned_bitsets_from_file(ireader.reader(inbed2))
		except:
			logging.error("invalid input: %s" % inbed2)
			sys.exit(1)
	else:
		logging.error("invalid input: %s" % inbed2)
		sys.exit(1)
	ret_lst = []
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
			ret_lst.append((chrom,start,end))
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


def bed_actual_size(*argv):
	'''
	Calculate the aggregated size of BED file.


	Parameters
	----------
	argv : list of genomic regions.
		Each argument can be a list, BED-like file, or a bigBed file. BED file
		can be regular, compressed, or remote file. The suffix of bigBed file
		must be one of ('.bb','.bigbed','.bigBed','.BigBed', '.BB',' BIGBED').

	Returns
	-------
	List of aggregated size.

	Examples
	--------
	>>> bed1 = [('chr1', 1, 10), ('chr1', 20, 35)]
	>>> bed2 = [('chr1',3, 15), ('chr1',20, 50)]
	>>> bed_actual_size(bed1, bed2)
	[24, 42]
	'''
	sizes = []
	for arg in argv:
		size = 0
		if type(arg) is list:
			for chrom,start,end in arg:
				size += (int(end) - int(start))
		elif type(arg) is str:
			for l in ireader.reader(arg):
				if l.startswith(('browser','#','track')):continue
				f = l.split()
				if len(f) < 3:
					logging.warning("invalid BED line: %s" % l)
					continue
				size += (int(f[2]) - int(f[1]))
		sizes.append(size)
	return sizes


def bed_counts(*argv):
	'''
	Calculate the number of regions in BED file.


	Parameters
	----------
	argv : list of genomic regions.
		Each argument can be a list, BED-like file, or a bigBed file. BED file
		can be regular, compressed, or remote file. The suffix of bigBed file
		must be one of ('.bb','.bigbed','.bigBed','.BigBed', '.BB',' BIGBED').

	Returns
	-------
	List of aggregated size.

	Examples
	--------
	>>> bed1 = [('chr1', 1, 10), ('chr1', 20, 35)]
	>>> bed2 = [('chr1',3, 15), ('chr1',20, 50), ('chr2',100,200)]
	>>> bed_counts(bed1, bed2)
	[2, 3]
	'''
	bed_counts = []
	for arg in argv:
		count = 0
		if type(arg) is list:
			count = len(arg)
		elif type(arg) is str:
			for l in ireader.reader(arg):
				if l.startswith(('browser','#','track')):continue
				f = l.split()
				if len(f) < 3:
					logging.warning("invalid BED line: %s" % l)
					continue
				count += 1
		bed_counts.append(count)
	return bed_counts


def bed_genomic_size(*argv):
	'''
	Calculate the *genomic size* of BED file (or list). Note, genomic_size <= actual_size.

	Parameters
	----------
	argv : list of genomic regions.
		Each argument can be a list, BED-like file, or a bigBed file. BED file
		can be regular, compressed, or remote file. The suffix of bigBed file
		must be one of ('.bb','.bigbed','.bigBed','.BigBed', '.BB',' BIGBED').

	Returns
	-------
	List of genomic size.

	Example
	-------
	>>> bed1 = [('chr1', 0, 100), ('chr1', 50, 150), ('chr1', 80, 180)]
	>>> bed_genomic_size(bed1)
	[180]
	>>> bed2 = [('chr1', 0, 100), ('chr2', 50, 150), ('chr3', 80, 180)]
	>>> bed_genomic_size(bed2)
	[300]
	>>> bed_genomic_size(bed1, bed2)
	[180, 300]
	'''
	union_sizes = []
	for arg in argv:
		if type(arg) is list:
			bitsets = binned_bitsets_from_list(arg )
		elif type(arg) is str:
			try:
				bitsets = binned_bitsets_from_file( ireader.reader(arg) )
			except:
				logging.error("Invalid input: %s" % arg)
				sys.exit(1)
		else:
			logging.error("Invalid input: %s" % arg)
			sys.exit(1)
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
			logging.error("invalid BED line: %s" % l)
			sys.exit(1)
		if (int(f[2]) - int(f[1])) < 0:
			logging.error("invalid BED line: %s" % l)
			sys.exit(1)
		regions.append(  (f[0], int(f[1]), int(f[2]))  )
	return regions


def compare_bed(inbed1, inbed2):

	"""
	Compare two BED files (or lists). This function is similar to Linux "comm" command.

	Parameters
	----------
	inbed1 : str or list
		Name of a BED file or list of BED regions, for example,
		[(chr1 100 200), (chr2 150  300), (chr2 1000 1200)]
	inbed2 : str or list
		Name of a BED file or list of BED regions, for example,
		[(chr1 100 200), (chr2 150  300), (chr2 1000 1200)]

	Returns
	-------
	bed1_uniq : list
		Genomic regions that are inbed1 unique (i.e., regions only present in inbed1 but
		do not overlap with any regions in inbed2).
	bed2_uniq : list
		Genomic regions that are inbed2 unique (i.e., regions only present in inbed2 but
		do not overlap with any regions in inbed1).
	common : list
		Genomic regions overlapped between inbed1 and inbed2. Note, the
		overlapped regions were merged. For example, (chr1 1 10) and (chr1 5 15)
		will be merged as (chr1 1 15).

	Note
	----
	Overlapped regions *within* input BED files (or lists) are merged before
	comparison.

	"""

	logging.info("Read and union BED file: \"%s\"" % inbed1)
	bed1_union = unionBed3(inbed1)
	#logging.info("Original regions of %s : %d" % (inbed1, len(inbed1)))
	logging.info("Unioned regions of \"%s\" : %d" % (inbed1, len(bed1_union)))

	logging.info("Read and union BED file: \"%s\"" % inbed2)
	bed2_union = unionBed3(inbed2)
	#logging.info("Original regions of %s : %d" % (inbed2, len(inbed2)))
	logging.info("Unioned regions of \"%s\" : %d" % (inbed2, len(bed2_union)))

	logging.info("Merge BED files \"%s\" and \"%s\"" % (inbed1, inbed2))
	bed12_union = unionBed3(bed1_union + bed2_union)
	logging.info("Unioned regions of two BED files : %d" % len(bed12_union))


	logging.info("Build interval tree for unioned BED file: \"%s\"" % inbed1)
	maps1 = {}
	for (ichr1, istart1, iend1) in bed1_union:
		if ichr1 not in maps1:
			maps1[ichr1] = Intersecter()
		maps1[ichr1].add_interval( Interval(istart1, iend1))

	logging.info("Build interval tree for unioned BED file: \"%s\"" % inbed2)
	maps2 = {}
	for (ichr2, istart2, iend2) in bed2_union:
		if ichr2 not in maps2:
			maps2[ichr2] = Intersecter()
		maps2[ichr2].add_interval( Interval(istart2, iend2))

	logging.info("Find common and specific regions ...")
	bed1_uniq = []
	bed2_uniq = []
	common = []
	for (chrom, start, end) in bed12_union:
		if chrom in maps1 and chrom in maps2:
			#found in maps1
			if len( maps1[chrom].find(start, end) ) > 0:
				#found in maps2
				if len( maps2[chrom].find(start, end) ) > 0:
					common.append((chrom, start, end))
				#not found in maps2
				else:
					bed1_uniq.append((chrom, start, end))
			#not found in maps1
			else:
				#found in maps2
				if len( maps2[chrom].find(start, end) ) > 0:
					bed2_uniq.append((chrom, start, end))
				#not found in maps2
				else:
					continue
		elif chrom in maps1:
			bed1_uniq.append((chrom, start, end))
		elif  chrom in maps2:
			bed2_uniq.append((chrom, start, end))
		else:
			continue
	logging.info("\"%s\" unique regions: %d" % (inbed1, len(bed1_uniq)))
	logging.info("\"%s\" unique regions: %d" % (inbed2, len(bed2_uniq)))
	logging.info("Common (overlapped) regions: %d" % len(common))
	return (bed1_uniq, bed2_uniq, common)


def compare_peak(inbed1, inbed2, na_label='NA'):

	"""
	Calculates peak-wise overlaps.

	Parameters
	----------
	inbed1 : str
		Name of a BED file.
	inbed2 : str
		Name of another BED file.
	na_label : str
		String label used to represent missing value.

	Returns
	-------
	None
	"""
	pattern = re.compile(".bed$", re.IGNORECASE)
	logging.info("Read and union BED file: \"%s\"" % inbed1)
	bed1_union = unionBed3(inbed1)
	logging.info("Unioned regions of \"%s\" : %d" % (inbed1, len(bed1_union)))

	logging.info("Read and union BED file: \"%s\"" % inbed2)
	bed2_union = unionBed3(inbed2)
	logging.info("Unioned regions of \"%s\" : %d" % (inbed2, len(bed2_union)))

	#logging.info("Merge BED files \"%s\" and \"%s\"" % (inbed1, inbed2))
	#bed12_union = unionBed3(bed1_union + bed2_union)
	#logging.info("Unioned regions of two BED files : %d" % len(bed12_union))


	logging.info("Build interval tree for unioned BED file: \"%s\"" % inbed1)
	maps1 = {}
	for (ichr1, istart1, iend1) in bed1_union:
		if ichr1 not in maps1:
			maps1[ichr1] = Intersecter()
		maps1[ichr1].add_interval( Interval(istart1, iend1))

	logging.info("Build interval tree for unioned BED file: \"%s\"" % inbed2)
	maps2 = {}
	for (ichr2, istart2, iend2) in bed2_union:
		if ichr2 not in maps2:
			maps2[ichr2] = Intersecter()
		maps2[ichr2].add_interval( Interval(istart2, iend2))

	#overlap bed file 1 with bed file 2
	logging.info("Calculate the overlap coefficient of each genomic region in %s ..." % inbed1)
	outfile_name1 = pattern.sub('_ovcoef.tsv', os.path.basename(inbed1))
	BED1OUT = open(outfile_name1, 'w')
	print('\t'.join(['chrom','start','end','ov_peaks_n','ov_bases_n', 'ov_bases_frac','ov_coef','ov_peaks_list']), file=BED1OUT)
	for chrom, start, end in bed1_union:
		try:
			bed_1_size = end - start
			if bed_1_size <= 0:
				logging.debug("Skip %s" % (chrom + ':' + str(start) + '-' + str(end)))
			bed_1_lst = [(chrom, start, end)]
			bed_2_size = 0
			bed_2_lst = []

			overlaps = maps2[chrom].find(start, end)
			#print (overlaps)
			if len(overlaps) == 0:
				print('\t'.join([str(i) for i in (chrom, start, end, 0, na_label,na_label,na_label, na_label)]), file=BED1OUT)
			else:
				for o in overlaps:
					bed_2_size += (o.end - o.start)
					bed_2_lst.append((chrom, o.start, o.end))
				overlap_size = bed_overlap_size(bed_1_lst, bed_2_lst)
				try:
					peak_ov_coef = overlap_size/(math.sqrt(bed_1_size * bed_2_size))
				except:
					peak_ov_coef = 0
				tmp = ','.join([i[0] + ':' + str(i[1]) + '-' + str(i[2]) for i in bed_2_lst])
				print('\t'.join([str(i) for i in (chrom, start, end, len(bed_2_lst), overlap_size, overlap_size/bed_1_size, peak_ov_coef, tmp)]), file=BED1OUT)
		except:
			print('\t'.join([str(i) for i in (chrom, start, end, len(bed_2_lst), na_label, na_label, na_label, na_label)]), file=BED1OUT)
	BED1OUT.close()

	#overlap bed file 2 with bed file 1
	logging.info("Calculate the overlap coefficient of each genomic region in %s ..." % inbed2)
	outfile_name2 = pattern.sub('_ovcoef.tsv', os.path.basename(inbed2))
	BED2OUT = open(outfile_name2, 'w')
	print('\t'.join(['chrom','start','end','ov_peaks_n','ov_bases_n', 'ov_bases_frac','ov_coef','ov_peaks_list']), file=BED2OUT)
	for chrom, start, end in bed2_union:
		try:
			bed_2_size = end - start
			if bed_2_size <= 0:
				logging.debug("Skip %s" % (chrom + ':' + str(start) + '-' + str(end)))
			bed_2_lst = [(chrom, start, end)]
			bed_1_size = 0
			bed_1_lst = []

			overlaps = maps1[chrom].find(start, end)
			if len(overlaps) == 0:
				print('\t'.join([str(i) for i in (chrom, start, end, 0, na_label,na_label,na_label, na_label)]), file=BED2OUT)
			else:
				for o in overlaps:
					bed_1_size += (o.end - o.start)
					bed_1_lst.append((chrom, o.start, o.end))
				overlap_size = bed_overlap_size(bed_2_lst, bed_1_lst)
				try:
					peak_ov_coef = overlap_size/(math.sqrt(bed_1_size * bed_2_size))
				except:
					peak_ov_coef = 0
				tmp = ','.join([i[0] + ':' + str(i[1]) + '-' + str(i[2]) for i in bed_1_lst])
				print('\t'.join([str(i) for i in (chrom, start, end, len(bed_2_lst), overlap_size, overlap_size/bed_2_size, peak_ov_coef, tmp)]), file=BED2OUT)
		except:
			print('\t'.join([str(i) for i in (chrom, start, end, len(bed_1_lst), na_label, na_label, na_label, na_label)]), file=BED2OUT)
	BED2OUT.close()


def cooccur_peak(inbed1, inbed2, inbed_bg, outfile, n_cut=1, p_cut=0.0):

	"""
	Evaluate if two peak sets are significantly oc-occurred or mutually exclusive.
	Using Fisher's exact test.

	Parameters
	----------
	inbed1 : str
		Name of a BED file.
	inbed2 : str
		Name of another BED file.
	inbed_bg : str
		Name of the background BED file (e.g., all promoters, all enhancers).
	outfile : str
		Name of the output file.
	n_cut : int, optional
		Threshold of overlap size. For example, the overlap size is 20 for these two
		regions ('chr1', 0, 100) and ('chr1', 80, 250).
		default = 1
	p_put : float, optional
		Threshold of overlap percentage. In the example above, the overlap percentage
		for ('chr1', 0, 100) is 20/100 = 0.2.
		default = 0.0
	Returns
	-------
	None
	"""

	logging.info("Read and union BED file: \"%s\"" % inbed1)
	bed1_union = unionBed3(inbed1)
	logging.info("Number of unioned regions : %d" % (len(bed1_union)))

	logging.info("Read and union BED file: \"%s\"" % inbed2)
	bed2_union = unionBed3(inbed2)
	logging.info("Number of unioned regions : %d" % (len(bed2_union)))

	if inbed_bg is None:
		logging.info("Merge two input BED files \"%s\" and \"%s\" as the background" % (inbed1, inbed2))
		background = unionBed3(bed1_union + bed2_union)
		logging.info("Number of unioned background regions : %d" % (len(background)))
	else:
		logging.info("Read and union BED file: \"%s\"" % inbed_bg)
		background = unionBed3(inbed_bg)
		logging.info("Number of unioned background regions  : %d" % len(background))


	logging.info("Build interval tree for : \"%s\"" % inbed1)
	maps1 = {}
	for (ichr1, istart1, iend1) in bed1_union:
		if ichr1 not in maps1:
			maps1[ichr1] = Intersecter()
		maps1[ichr1].add_interval( Interval(istart1, iend1))

	logging.info("Build interval tree for: \"%s\"" % inbed2)
	maps2 = {}
	for (ichr2, istart2, iend2) in bed2_union:
		if ichr2 not in maps2:
			maps2[ichr2] = Intersecter()
		maps2[ichr2].add_interval( Interval(istart2, iend2))

	#background regions will be divided into 4 categories
	bed1_only = 0
	bed2_only = 0
	cooccur = 0
	neither = 0
	OUT = open(outfile, 'w')

	results = {}
	for chrom, start, end in background:
		line = chrom + '\t' + str(start) + '\t' + str(end)
		bed1_flag = False
		bed2_flag = False

		if (chrom not in maps1) and (chrom not in maps2):
			pass
		elif (chrom not in maps1) and (chrom in maps2):
			bed2_overlaps = maps2[chrom].find(start, end)
			if len(bed2_overlaps) == 0:
				pass
			else:
				bed2_overlap_lst = []
				for o in bed2_overlaps:
					bed2_overlap_lst.append((chrom, o.start, o.end))
				bed2_overlap_size = bed_overlap_size([(chrom, start, end)], bed2_overlap_lst)
				bed2_genomic_size = bed_genomic_size(bed2_overlap_lst)
				try:
					bed2_overlap_ratio = bed2_overlap_size/bed2_genomic_size
				except:
					bed2_overlap_ratio = 0
				if bed2_overlap_size >= n_cut and bed2_overlap_ratio >= p_cut:
					bed2_flag = True
		elif (chrom in maps1) and (chrom not in maps2):
			bed1_overlaps = maps1[chrom].find(start, end)
			if len(bed2_overlaps) == 0:
				pass
			else:
				bed1_overlap_lst = []
				for o in bed1_overlaps:
					bed1_overlap_lst.append((chrom, o.start, o.end))
				bed1_overlap_size = bed_overlap_size([(chrom, start, end)], bed1_overlap_lst)
				bed1_genomic_size = bed_genomic_size(bed1_overlap_lst)
				try:
					bed1_overlap_ratio = bed1_overlap_size/bed1_genomic_size
				except:
					bed1_overlap_ratio = 0
				if bed1_overlap_size >= n_cut and bed1_overlap_ratio >= p_cut:
					bed1_flag = True
		else:
			#overlaps with inbed1
			bed1_overlaps = maps1[chrom].find(start, end)
			#overlaps with inbed2
			bed2_overlaps = maps2[chrom].find(start, end)
			if len(bed1_overlaps) == 0:
				if len(bed2_overlaps) == 0:
					pass
				elif len(bed2_overlaps) > 0:
					bed2_overlap_lst = []
					for o in bed2_overlaps:
						bed2_overlap_lst.append((chrom, o.start, o.end))
					bed2_overlap_size = bed_overlap_size([(chrom, start, end)], bed2_overlap_lst)
					bed2_genomic_size = bed_genomic_size(bed2_overlap_lst)
					try:
						bed2_overlap_ratio = bed2_overlap_size/bed2_genomic_size
					except:
						bed2_overlap_ratio = 0
					if bed2_overlap_size >= n_cut and bed2_overlap_ratio >= p_cut:
						bed2_flag = True
			elif len(bed1_overlaps) > 0:
				bed1_overlap_lst = []
				for o in bed1_overlaps:
					bed1_overlap_lst.append((chrom, o.start, o.end))
				bed1_overlap_size = bed_overlap_size([(chrom, start, end)], bed1_overlap_lst)
				bed1_genomic_size = bed_genomic_size(bed1_overlap_lst)
				try:
					bed1_overlap_ratio = bed1_overlap_size/bed1_genomic_size
				except:
					bed1_overlap_ratio = 0
				if bed1_overlap_size >= n_cut and bed1_overlap_ratio >= p_cut:
					bed1_flag = True

				if len(bed2_overlaps) == 0:
					pass
				elif len(bed2_overlaps) > 0:
					bed2_overlap_lst = []
					for o in bed2_overlaps:
						bed2_overlap_lst.append((chrom, o.start, o.end))
					bed2_overlap_size = bed_overlap_size([(chrom, start, end)], bed2_overlap_lst)
					bed2_genomic_size = bed_genomic_size(bed2_overlap_lst)
					try:
						bed2_overlap_ratio = bed2_overlap_size/bed2_genomic_size
					except:
						bed2_overlap_ratio = 0
					if bed2_overlap_size >= n_cut and bed2_overlap_ratio >= p_cut:
						bed2_flag = True
		if bed1_flag:
			if bed2_flag:
				cooccur += 1
				print (line + '\tCooccur', file=OUT)
			else:
				bed1_only += 1
				print (line + '\t%s_only' % os.path.basename(inbed1), file=OUT)
		else:
			if bed2_flag:
				bed2_only += 1
				print (line + '\t%s_only' % os.path.basename(inbed2), file=OUT)
			else:
				neither += 1
				print (line + '\tNeither', file=OUT)
	#print (bed1_only, bed2_only, cooccur, neither)
	results['bed1_name'] = os.path.basename(inbed1)
	results['bed2_name'] = os.path.basename(inbed2)
	results['bed1+,bed2-'] = bed1_only
	results['bed1-,bed2+'] = bed2_only
	results['bed1+,bed2+'] = cooccur
	results['bed1-,bed2-'] = neither

	results['Jaccard index'] = cooccur/(cooccur + neither + bed1_only + bed2_only)
	if bed1_only > bed2_only:
		table = np.array([[neither, bed1_only], [bed2_only, cooccur]])
	else:
		table = np.array([[neither, bed2_only], [bed1_only, cooccur]])
	print (table)
	oddsr,p = fisher_exact(table, alternative='greater')
	results['odds-ratio'] = oddsr
	results['p-value'] = p
	return pd.Series(data=results, index=['bed1_name', 'bed2_name','bed1+,bed2-', 'bed1-,bed2+', 'bed1+,bed2+', 'bed1-,bed2-','odds-ratio', 'p-value', 'Jaccard index'], name = "Fisher's exact test result")


def bedtofile(bed_list, bed_file):
	''' Save list of genomic regions to file'''
	OUT = open(bed_file,'w')
	for tmp in bed_list:
		print ('\t'.join([str(i) for i in tmp]), file=OUT)
	OUT.close()

def is_overlap(chr1, st1, end1, chr2, st2, end2):
	'''
	Check if two regios are overlap.

	Parameters
	----------
	chr1 : str
		Chromosome ID of the first genomic region
	st1 : int
		Start coordinate of the first genomic region
	end1 : int
		End coordinate of the first genomic region
	chr2 : str
		Chromosome ID of the second genomic region
	st2 : int
		Start coordinate of the second genomic region
	end2 : int
		End coordinate of the second genomic region
	Return
	------
		0: non overlap
		positive integer [1,): overlap
	'''
	#genome coordinate is left-open, right-close.
	st1 = st1 +1
	end1 = end1
	st2 = st2 +1
	end2 = end2
	if chr1 != chr2:
		return 0
	else:
		return len(range(max(st1, st2), min(end1, end2)+1))

def srogcode(lst1, lst2):
	"""
	Determine the spatial relations of genomic regions (SROG)

	Parameters
	----------
	lst1 : tuple
		A tuple of genomic interval. Eg. ('chr1', 1, 100).
	lst2 : tuple
		A tuple of genomic interval. Eg. ('chr1', 15, 60)..

	Returns
	-------
	return_code : str
		Return one of ('disjoint','touch','equal','overlap','contain','within').

	"""
	return_code = ''
	try:
		chrom_1 = lst1[0]
		start_1 = int(lst1[1])
		end_1 = int(lst1[2])
		chrom_2 = lst2[0]
		start_2 = int(lst2[1])
		end_2 = int(lst2[2])
	except:
		return_code = 'unknown'
		return return_code
	ov_size = is_overlap(chrom_1, start_1, end_1, chrom_2, start_2, end_2)
	if ov_size == 0:
		if chrom_1 != chrom_2:
			return_code = 'disjoint'
		else:
			if start_1 == end_2 or end_1 == start_2:
				return_code = 'touch'
			else:
				return_code = 'disjoint'
	elif ov_size > 0:
		if start_1 == start_2 and end_1 == end_2:
			return_code = 'equal'
		elif start_1 >= start_2 and end_1 < end_2:
			return_code = 'within'
		elif start_1 > start_2 and end_1 <= end_2:
			return_code = 'within'
		elif start_1 > start_2 and end_1 < end_2:
			return_code = 'within'
		elif start_1 <= start_2 and end_1 > end_2:
			return_code = 'contain'
		elif start_1 < start_2 and end_1 >= end_2:
			return_code = 'contain'
		elif start_1 < start_2 and end_1 > end_2:
			return_code = 'contain'
		else:
			return_code = 'overlap'
	return return_code


def srog_peak(inbed1, inbed2, outfile, n_up = 1, n_down = 1):

	"""
	Calculates SROG code for each region in inbed1

	Parameters
	----------
	inbed1 : str
		Name of a BED file.
	inbed2 : str
		Name of another BED file.
	outfile : str
		Name of output file.

	Returns
	-------
	None
	"""

	maps = {}
	OUTPUT = open(outfile, 'w')
	srog_summary = {'disjoint':0, 'overlap':0, 'contain':0, 'within':0, 'touch':0, 'equal':0, 'other':0}
	logging.info("Build interval tree from file: \"%s\"" % inbed2)
	for l in ireader.reader(inbed2):
		if l.startswith(('browser','#','track')):continue
		f = l.split()
		if len(f) < 3:
			logging.warning("Invalid BED line (Requires at least 3 columns: chrom, start, end): %s" % l)
			continue
		chrom, start, end = f[0], int(f[1]), int(f[2])
		if start > end:
			logging.warning("invalid BED line (start > end): %s" % l)
			continue

		#try to get name from the 4th column
		try:
			name = f[3]
		except:
			name = f[0] + ':' + f[1] + '-' + f[2]

		#try to get strand from the 6th column
		try:
			strandness = f[5]
			if strandness not in ('+','-'):
				strandness = '+'
		except:
			strandness = '+'

		if chrom not in maps:
			maps[chrom] = Intersecter()
		maps[chrom].add_interval( Interval(start, end, value = name, strand = strandness))

	logging.info("Reading BED file: \"%s\"" % inbed1)
	for l in ireader.reader(inbed1):
		if l.startswith(('browser','#','track')):continue
		f = l.split()
		if len(f) < 3:
			logging.warning("Invalid BED line (Requires at least 3 columns: chrom, start, end): %s" % l)
			continue
		chrom, start, end = f[0], int(f[1]), int(f[2])
		if start > end:
			logging.warning("Invalid BED line (start > end): %s" % l)
			continue

		#try to get name from the 4th column
		try:
			name = f[3]
		except:
			name = f[0] + ':' + f[1] + '-' + f[2]

		#try to get strand from the 6th column
		try:
			strandness = f[5]
			if strandness not in ('+','-'):
				strandness = '+'
		except:
			strandness = '+'

		if chrom not in maps:
			srog_summary['disjoint'] += 1
			print (l + '\t' + 'NA' + '\t' + 'NA', file=OUTPUT)
			continue

		overlaps = maps[chrom].find(start, end)
		if len(overlaps) == 0:
			srog_summary['disjoint'] += 1
			up_interval =  maps[chrom].upstream_of_interval(Interval(start, end, strand = strandness), num_intervals = n_up)
			down_interval =  maps[chrom].downstream_of_interval(Interval(start, end, strand = strandness), num_intervals = n_down)
			if len(up_interval) == 0:
				up_interval_name = 'NA'
			else:
				up_interval_name = str(up_interval[0].value)
			if len(down_interval) == 0:
				down_interval_name = 'NA'
			else:
				down_interval_name = str(down_interval[0].value)
			print (l + '\t' + 'disjoint' + '\t' + 'UpInterval=' + up_interval_name + ',' + 'DownInterval=' + down_interval_name, file=OUTPUT)
		else:
			srog_codes = []
			target_names = []
			for o in overlaps:
				tmp = srogcode((chrom, start, end), (chrom, o.start, o.end))
				srog_codes.append(tmp)
				target_names.append(o.value)
			print (l + '\t' + ','.join(srog_codes) + '\t' + ','.join(target_names), file=OUTPUT)
			for code in srog_codes:
				srog_summary[code] += 1
	return pd.Series(data=srog_summary)


if __name__=='__main__':
	#(a, b, common) = compare_bed(sys.argv[1], sys.argv[2])
	#bedtofile(a,'a')
	#bedtofile(b,'b')
	#bedtofile(common,'common')

	a = cooccur_peak(sys.argv[1], sys.argv[2], sys.argv[3])
	print (a)