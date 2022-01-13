#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 13:06:22 2022

@author: m102324
"""

import sys
from cobindability.BED import bed_actual_size, bed_genomic_size, bed_overlap_size, bedtolist, bed_counts
from cobindability.coefcal import ov_coef, ov_jaccard, ov_ss, ov_sd
from os.path import basename
import logging
import numpy as np
import pandas as pd
from random import sample
from cobindability import version


__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "MIT"
__version__ = version.version
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Development"

def jindex(union, overlap):
	"""
	Calculate Jaccard's index (a.k.a, Jaccard similarity coefficient).

	Parameters
	----------

	union : int
		Number of total bases (non-redundant) in A and B.
	overlap : int
		Number of bases shared (i.e., overlapped) between A and B.

	Returns
	-------
	jaccard_index : FLOAT
		The observed Jaccard's index or Jaccard's coefficient.
	"""

	jaccard_index = 0
	if overlap < 0 or union < 0:
		logging.error("Invalid parameters!")
		return None
	if overlap > union:
		logging.error("Intersection larger than union!")
		return None
	if union == 0:
		logging.warning("Cannot calculate the observed Jaccard index (coefficient)!")
		return None
	jaccard_index = overlap/union
	if jaccard_index >= 0 and jaccard_index <= 1:
		return jaccard_index
	else:
		return None

def jindex_exp(size_a, size_b, size_bg = 1.4E9):

	"""
	Calculate *expected* Jaccard's index (a.k.a, Jaccard similarity coefficient).

	Parameters
	----------
	size_a : int
		Number of total bases (non-redundant) in interval set A.
	size_b : int
		Number of total bases (non-redundant) in interval set B.
	size_bg : int, optional
		The effective background genome size.
		default = 1.4E9

	Returns
	-------
	jaccard_index_exp : FLOAT
		The expected Jaccard's index or Jaccard's coefficient.
	"""
	jaccard_index_exp = 0
	if min(size_a, size_b, size_bg) < 0 or size_bg == 0:
		logging.error("Invalid parameters!")
		return None

	if max(size_a, size_b, size_bg) != size_bg:
		logging.error("Invalid parameters!")
		return None

	expected_overlaps = size_a * size_b/size_bg
	union_size = size_a + size_b - expected_overlaps
	if union_size <= 0:
		return None
	jaccard_index_exp = expected_overlaps/union_size
	if jaccard_index_exp >= 0 and jaccard_index_exp <= 1:
		return jaccard_index_exp
	else:
		return None


def sdindex(size_a, size_b, size_overlap):

	"""
	Calculate the Sørensen–Dice index (a.k.a, Sørensen–Dice coefficient).

	Parameters
	----------
	size_a : int
		Number of total bases (non-redundant) in interval set A.
	size_b : int
		Number of total bases (non-redundant) in interval set B.
	size_overlap : int
		Number of overlapped bases between A and B.

	Returns
	-------
	sd_index : FLOAT
		The expected Jaccard's index or Jaccard's coefficient.
	"""
	sd_index = 0
	if min(size_a, size_b, size_overlap) < 0:
		logging.error("Invalid parameters!")
		return None
	if size_a + size_b == 0:
		logging.error("Invalid parameters!")
		return None
	if size_overlap > size_a or size_overlap > size_b:
		logging.error("Invalid parameters!")
		return None

	sd_index = 2.0*size_overlap/(size_a + size_b)
	if sd_index >= 0 and sd_index <= 1:
		return sd_index
	else:
		return None


def sdindex_exp(size_a, size_b, size_bg = 1.4E9):

	"""
	Calculate *expected* Sørensen–Dice index (a.k.a, Sørensen–Dice coefficient).

	Parameters
	----------
	size_a : int
		Number of total bases (non-redundant) in interval set A.
	size_b : int
		Number of total bases (non-redundant) in interval set B.
	size_bg : int, optional
		The effective background genome size.
		default = 1.4E9

	Returns
	-------
	sd_index_exp : FLOAT
		The expected Sørensen–Dice index or Sørensen–Dice coefficient.
	"""
	sd_index_exp = 0
	if min(size_a, size_b, size_bg) < 0 or size_bg == 0:
		logging.error("Invalid parameters!")
		return None

	if max(size_a, size_b, size_bg) != size_bg:
		logging.error("Invalid parameters!")
		return None

	expected_overlaps = size_a * size_b/size_bg

	sd_index_exp = 2.0*expected_overlaps/(size_a + size_b)
	if sd_index_exp >= 0 and sd_index_exp <= 1:
		return sd_index_exp
	else:
		return None




def cal_jaccard_coef(file1, file2, n_draws=10, fraction = 0.5, bg_size = 1.4e9):
	"""
	Calculate Jaccard's index (a.k.a, Jaccard similarity coefficient) between two BED files.

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
	n_draws : int
		Times of drawing samples with replacement. Set to '0' to turn off bootstraping.
	fraction : float, optional
		The fraction of subsample. The default is 0.85 (85% of the orignal genomic reginos will be selected).
	bg_size : int, optional
		The effective background genome size. About 1.4Gb of the human genome are. The default is 1.4e9.

	Returns
	-------
	pandas.Series
		Pandas Series of 'coef_obs', 'coef_exp','coef_ratio', 'coef_ratio_low', 'coef_ratio_high'.
		'coef_obs' : Observed overlap coefficient between two BED files.
		'coef_exp' : Expected overlap coefficient between two BED files.
		'coef_ratio' : Ratio between 'coef_obs' and 'coef_exp'.
		'coef_ratio_low' : The lower bound of 95% confidence interval of 'coef_ratio'.
		'coef_ratio_high' : The upper bound of 95% confidence interval of 'coef_ratio'.

	"""

	results = {}

	file1_lst = bedtolist (file1)
	file2_lst = bedtolist (file2)

	file1_name = basename(file1)
	file2_name = basename(file2)

	#calculate interval counts
	logging.debug("Calculating bed counts ...")
	totalCount1, totalCount2 = bed_counts(file1, file2)
	results[file1_name + '.count'] = totalCount1
	results[file2_name + '.count'] = totalCount2


	#calculate overall overlap coef
	logging.info("Calculating the overall Jaccard's coefficient ...")
	(uniqBase1,uniqBase2) = bed_genomic_size(file1_lst, file2_lst)


	overlapBases = bed_overlap_size(file1_lst, file2_lst)
	overlapBases_exp = uniqBase1*uniqBase2/bg_size
	[unionBases] = bed_genomic_size(file1_lst + file2_lst)

	results[file1_name + '.uniq_size'] = uniqBase1
	results[file2_name + '.uniq_size'] = uniqBase2
	results['union_size'] = unionBases
	results['overlap_size'] = overlapBases
	results['Jaccard_coef'] = ov_jaccard(uniqBase1, uniqBase2, overlapBases)
	results['Jaccard_coef_exp'] = ov_jaccard(uniqBase1, uniqBase2, overlapBases_exp)
	results['Jaccard_distance'] = 1 - results['Jaccard_coef']


	#logging.info("Calculating the overall Sørensen–Dice coefficient ...")
	#results['Sørensen–Dice_coefficient(obs)'] = sdindex(uniqBase1, uniqBase2, overlapBases)
	#results['Sørensen–Dice_coefficient(exp)'] = sdindex_exp(uniqBase1, uniqBase2, size_bg = bg_size)

	if n_draws > 0:
		if fraction > 0 and fraction < 1:
			resample_size1 = int(totalCount1 * fraction)
			resample_size2 = int(totalCount2 * fraction)
		else:
			logging.error("Fraction must be > 0 and < 1.")
			sys.exit(0)
		logging.debug("Bootstraping is on. Iterate %d times. " % n_draws)
		tmp = []
		for i in range(n_draws):
			logging.info("Bootstraping iteration %d ..." % i)
			sample_1 = sample(file1_lst, resample_size1)
			sample_2 = sample(file2_lst, resample_size2)
			overlapBases = bed_overlap_size(sample_1, sample_2)
			overlapBases = overlapBases * (1/fraction)
			(sample1_size, sample2_size) = bed_genomic_size(sample_1, sample_2)
			sample_ovcoef = ov_jaccard(sample1_size, sample2_size, overlapBases)

			tmp.append(sample_ovcoef)


		ci_lower = np.percentile(np.array(tmp), 2.5)
		ci_upper = np.percentile(np.array(tmp), 97.5)
		results['Jaccard_coef_obs_low'] = ci_lower
		results['Jaccard_coef_obs_high'] = ci_upper
	else:
		logging.info("Bootstraping is off ...")
		results['Jaccard_coef_obs_low'] = np.nan
		results['Jaccard_coef_obs_high'] = np.nan

	return pd.Series(data=results)

if __name__ == '__main__':

	logging.basicConfig(format = "%(asctime)s [%(levelname)s]  %(message)s",datefmt='%Y-%m-%d %I:%M:%S', level=logging.INFO)

	a = cal_jaccard_coef(sys.argv[1], sys.argv[2], n_draws=0, fraction = 0.8)
	print (a)
