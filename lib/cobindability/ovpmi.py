#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 13:11:02 2022

@author: m102324
"""

import sys
#import math
from cobindability.BED import bed_genomic_size, bed_overlap_size, bedtolist
import logging
import numpy as np
import pandas as pd
from cobindability import version

__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "MIT"
__version__ = version.version
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Development"


def pmi_value(px, py, pxy):
	if px < 0 or px >1:
		logging.error("Probability must be in (0,1)")
	if py < 0 or py >1:
		logging.error("Probability must be in (0,1)")
	if pxy < 0 or pxy >1:
		logging.error("Probability must be in (0,1)")
	if pxy == 0:
		return -np.inf
	else:
		return np.log2(pxy) - np.log2(px) - np.log2(py)

def npmi_value(px, py, pxy):
	if px < 0 or px >1:
		logging.error("Probability must be in (0,1)")
	if py < 0 or py >1:
		logging.error("Probability must be in (0,1)")
	if pxy < 0 or pxy >1:
		logging.error("Probability must be in (0,1)")

	if pxy == 0:
		return -1
	else:
		return np.log2(px*py)/np.log2(pxy) - 1


def cal_pmi(file1, file2, bg_size = 1424655930):
	"""

	Parameters
	----------
	file1 : str
		Genomic regions in BED (Browser Extensible Data, https://genome.ucsc.edu/FAQ/FAQformat.html#format1), BED-like or BigBed format.
		The BED-like format includes 'bed3','bed4','bed6','bed12','bedgraph','narrowpeak', 'broadpeak','gappedpeak'. BED and BED-like
		format can be plain text, compressed (.gz, .z, .bz, .bz2, .bzip2) or remote (http://, https://, ftp://) files. Do not compress
		BigBed foramt. BigBed file can also be a remote file.
	file2 : str
		Genomic regions in BED (Browser Extensible Data, https://genome.ucsc.edu/FAQ/FAQformat.html#format1), BED-like or BigBed format.
		The BED-like format includes 'bed3','bed4','bed6','bed12','bedgraph','narrowpeak', 'broadpeak','gappedpeak'. BED and BED-like
		format can be plain text, compressed (.gz, .z, .bz, .bz2, .bzip2) or remote (http://, https://, ftp://) files. Do not compress
		BigBed foramt. BigBed file can also be a remote file.
	bg_size : int, optional
		The effective background genome size. About 1.4Gb of the human genome are. The default is 1424655930.

	Returns
	-------
	List
		(PMI and NMPI)


	"""
	results = {}
	file1_lst = bedtolist (file1)
	file2_lst = bedtolist (file2)


	#calculate bed's genomic size
	(uniqBase1,uniqBase2) = bed_genomic_size(file1_lst, file2_lst)

	#calculate overall overlap coef
	overlapBases = bed_overlap_size(file1_lst, file2_lst)

	results['A_name'] = file1
	results['B_name'] = file2
	results['A'] = uniqBase1/bg_size
	results['not_A'] = (bg_size - uniqBase1)/bg_size
	results['B'] = uniqBase2/bg_size
	results['not_B'] = (bg_size - uniqBase2)/bg_size
	results['A_not_B'] = (uniqBase1 - overlapBases)/ bg_size
	results['B_not_A'] = (uniqBase2 - overlapBases)/bg_size
	results['A_and_B'] = overlapBases/bg_size
	results['Neither_A_nor_B'] = (bg_size - uniqBase1 - uniqBase2 + overlapBases)/bg_size

	results['A_and_B.pmi'] = pmi_value(results['A'], results['B'], results['A_and_B'])
	results['A_and_B.npmi'] = npmi_value(results['A'], results['B'], results['A_and_B'])

	results['A_not_B.pmi'] = pmi_value(results['A'], results['not_B'], results['A_not_B'])
	results['A_not_B.npmi'] = npmi_value(results['A'], results['not_B'], results['A_not_B'])

	results['B_not_A.pmi'] = pmi_value(results['B'], results['not_A'], results['B_not_A'])
	results['B_not_A.npmi'] = npmi_value(results['B'], results['not_A'], results['B_not_A'])

	return pd.Series(data=results,index=['A_name','B_name','A','not_A','B','not_B','A_not_B','B_not_A','A_and_B','Neither_A_nor_B','A_and_B.pmi','A_and_B.npmi', 'A_not_B.pmi', 'A_not_B.npmi', 'B_not_A.pmi', 'B_not_A.npmi'])


if __name__ == '__main__':

	b = cal_pmi(sys.argv[1], sys.argv[2])
	print (b)
