#!/usr/bin/env python

import sys
import math
from cobindability.BED import bed_actual_size, bed_genomic_size, bed_overlap_size, bedtolist
import logging
import numpy as np
import pandas as pd
from random import sample

__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="0.0.1"
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



def cal_overlap_coef(file1, file2, n_draws, size = 0.85, bg_size = 1424655930):
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
	n_draws : int
		Times of drawing samples with replacement. Set to '0' to turn off bootstraping.
	size : float, optional
		The size of subsample. The default is 0.85 (85% of the orignal genomic reginos will be selected).
	bg_size : int, optional
		The effective background genome size. About 1.4Gb of the human genome are. The default is 1424655930.

	Returns
	-------
	pandas.Series
		DESCRIPTION.

	"""

	results = {}

	file1_lst = bedtolist (file1)
	file2_lst = bedtolist (file2)

	#calculate bed total size and bed count
	totalBase1, totalBase2, totalCount1, totalCount2 = bed_actual_size(file1, file2)

	#calculate bed's genomic size
	(uniqBase1,uniqBase2) = bed_genomic_size(file1_lst, file2_lst)


	#calculate overall overlap coef
	overlapBases = bed_overlap_size(file1_lst, file2_lst)

	logging.debug("Calculating the overall observed overlapping coefficient ...")
	overlapCoef_obs = obs_coef(uniqBase1, uniqBase2, overlapBases)

	logging.debug("Calculating the overall expected overlapping coefficient ...")
	overlapCoef_exp = exp_coef(uniqBase1, uniqBase2, bg_size = bg_size)
	try:
		ratio = overlapCoef_obs/overlapCoef_exp
	except:
		logging.warning("Cannot calculate the ratio between observed coefficient and expected coefficient, set to 0.0!")
		ratio = 0.0
	results['coef_obs'] = overlapCoef_obs
	results['coef_exp'] = overlapCoef_exp
	results['coef_ratio'] = ratio
	if n_draws > 0:
		#bootstraping
		logging.debug("Bootstraping is on. Iterate %d times. " % n_draws)
		tmp = []
		for i in range(n_draws):
			logging.debug("Bootstraping iteration: %d ..." % i)
			sample_1 = sample(file1_lst, int(totalCount1*size))
			sample_2 = sample(file2_lst, int(totalCount2*size))
			(uniqBase1,uniqBase2) = bed_genomic_size(sample_1, sample_2)
			overlapBases = bed_overlap_size(sample_1, sample_2)

			overlapCoef_obs = obs_coef(uniqBase1, uniqBase2, overlapBases)
			overlapCoef_exp = exp_coef(uniqBase1, uniqBase2, bg_size = bg_size)

			#bs_coefs.append(overlapCoef)
			try:
				tmp.append(overlapCoef_obs/overlapCoef_exp)
			except:
				logging.warning("Cannot calculate the ratio between observed coefficient and expected coefficient, set to 0.0!")
				tmp.append(0,0)

		ci_lower = np.percentile(np.array(tmp), 2.5)
		ci_upper = np.percentile(np.array(tmp), 97.5)
		results['coef_ratio_low'] = ci_lower
		results['coef_ratio_high'] = ci_upper
	else:
		logging.info("Bootstraping is off ...")
		results['coef_ratio_low'] = np.nan
		results['coef_ratio_high'] = np.nan

	return pd.Series(data=results, index=['coef_obs', 'coef_exp','coef_ratio', 'coef_ratio_low', 'coef_ratio_high'])

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
	#return np.log2(pxy) - np.log2(px) - np.log2(py)


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
	#a = cal_overlap_coef(sys.argv[1], sys.argv[2], n_draws=0)
	#print (a)

	b = cal_pmi(sys.argv[1], sys.argv[2])
	print (b)
