#!/usr/bin/env python
# encoding: utf-8

"""
read bigBed files.

"""


import pyBigWig


def bbreader(fname):
	"""
	read bigBed file
	"""
	bb = pyBigWig.open(fname)
	chrom_dict = bb.chroms()
	for chr in chrom_dict:
		for start, end, score in bb.entries(chr, 0, chrom_dict[chr]):
			yield(chr + '\t' + str(start) + '\t' + str(end) + '\t' + score)

#test
"""
import sys
if __name__=='__main__':
	a = bbreader(sys.argv[1])
	for i in a:
		print (i)
"""