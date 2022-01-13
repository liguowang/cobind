#!/usr/bin/env python

import sys
import math
from cobindability.BED import bed_actual_size, bed_genomic_size, bed_overlap_size, bedtolist, bed_counts
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


def obs_coef(set1, set2, overlap):
	"""
	Calculate the *observed* overlapping coefficient between set1 and set2.

	Parameters
	----------
	set1 : int
		Number of nucleotides in genome region set1.
	set2 : int
		Number of nucleotides in genome region set2.
	overlap : int
		Number of nucleotides shared (i.e., overlapped) between set1 and set2.

	Returns
	-------
	observed_coef : FLOAT
		Overlapping coefficient.
	"""

	observed_coef = 0
	try:
		observed_coef = overlap/(math.sqrt(set1 * set2))
	except:
		logging.warning("Cannot calculate the observed overlapping coefficient!")
		pass
	return observed_coef

def exp_coef(set1, set2, bg_size = 1424655930):

	"""
	Calculate the *expected* overlapping coefficient between set1 and set2.

	Parameters
	----------
	set1 : int
		Number of nucleotides in genome region set1.
	set2 : int
		Number of nucleotides in genome region set2.
	bg_size : int, optional
		The effective background genome size. About 1.4Gb of the human genome are
		bound by transcription factors.

	Returns
	-------
	observed_coef : float
		The expected overlapping coefficient.
	"""
	expected_coef = 0
	expected_overlap = set1 * set2/bg_size
	try:
		expected_coef = expected_overlap/(math.sqrt(set1 * set2))
	except:
		logging.warning("Cannot calculate the expected overlapping coefficient!")
		pass
	return expected_coef



def cal_overlap_coef(file1, file2, n_draws, fraction = 0.75, bg_size = 1424655930):
	"""
	Calculate the *overall* overlap coefficient between two BED files.

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
		The effective background genome size. About 1.4Gb of the human genome are. The default is 1424655930.

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

	#calculate bed count
	logging.debug("Calculating bed counts ...")
	totalCount1, totalCount2 = bed_counts(file1, file2)
	results[file1_name + '.count'] = totalCount1
	results[file2_name + '.count'] = totalCount2

	#calculate bed total size
	logging.debug("Calculating total (aggregated) size ...")
	totalBase1, totalBase2 = bed_actual_size(file1, file2)
	results[file1_name + '.totalSize'] = totalBase1
	results[file2_name + '.totalSize'] = totalBase2

	#calculate bed's genomic (unique) size
	logging.debug("Calculating unique size ...")
	(uniqBase1,uniqBase2) = bed_genomic_size(file1_lst, file2_lst)
	results[file1_name + '.uniqSize'] = uniqBase1
	results[file2_name + '.uniqSize'] = uniqBase2

	#calculate overall overlap coef
	logging.debug("Calculating overlap size ...")
	overlapBases = bed_overlap_size(file1_lst, file2_lst)
	results['overlapSize'] = overlapBases

	logging.debug("Calculating the observed overall overlapping coefficient ...")
	overlapCoef_obs = obs_coef(uniqBase1, uniqBase2, overlapBases)
	#print (uniqBase1, uniqBase2, overlapBases)

	logging.debug("Calculating the expected overall overlapping coefficient ...")
	overlapCoef_exp = exp_coef(uniqBase1, uniqBase2, bg_size = bg_size)

	results['coef.overlap'] = overlapCoef_obs
	results['coef.overlap_exp'] = overlapCoef_exp
	if n_draws > 0:
		#bootstraping
		logging.debug("Bootstraping is on. Iterate %d times. " % n_draws)
		tmp = []
		for i in range(n_draws):
			logging.debug("Bootstraping iteration: %d ..." % i)
			sample_1 = sample(file1_lst, int(totalCount1*fraction))
			sample_2 = sample(file2_lst, int(totalCount2*fraction))
			(uniqBase1,uniqBase2) = bed_genomic_size(sample_1, sample_2)
			overlapBases = bed_overlap_size(sample_1, sample_2)

			overlapCoef_obs = obs_coef(uniqBase1, uniqBase2, overlapBases* (1/fraction))
			tmp.append(overlapCoef_obs)
		ci_lower = np.percentile(np.array(tmp), 2.5)
		ci_upper = np.percentile(np.array(tmp), 97.5)
		results['coef.overlap_low'] = ci_lower
		results['coef.overlap_high'] = ci_upper
	else:
		logging.info("Bootstraping is off ...")
		results['coef_ratio_low'] = np.nan
		results['coef_ratio_high'] = np.nan

	return pd.Series(data=results)



if __name__ == '__main__':
	a = cal_overlap_coef(sys.argv[1], sys.argv[2], n_draws=10)
	print (a)

