#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 14:17:50 2022

@author: m102324
"""

import sys
import logging
import pandas as pd
from cobindability.BED import bed_overlap_size, bed_to_list, bed_info
from cobindability.coefcal import ov_coef, ov_jaccard, ov_ss, ov_sd, pmi_value, npmi_value
from cobindability import version


__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "MIT"
__version__ = version.version
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Development"


def ov_stats(file1, file2, bg_size = 1400000000):
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
	draws_n : int
		Times of drawing samples with replacement. Set to '0' to turn off bootstraping.
	draw_frac : float, optional
		The size of subsample. The default is 0.85 (85% of the orignal genomic reginos will be selected).
	bg_size : int, optional
		The effective background genome size. About 1.4Gb of the human genome are. The default is 1424655930.

	Returns
	-------
	pandas.Series
		Pandas Series of 'coef_obs', 'coef_exp','coef_ratio', 'coef_ratio_low', 'coef_ratio_high'.
		'coef_obs' : Observed collocation coefficient between two BED files.
		'coef_exp' : Expected collocation coefficient between two BED files.
		'coef_ratio' : Ratio between 'coef_obs' and 'coef_exp'.
		'coef_ratio_low' : The lower bound of 95% confidence interval of 'coef_ratio'.
		'coef_ratio_high' : The upper bound of 95% confidence interval of 'coef_ratio'.


	"""
	results = {}
	file1_lst = bed_to_list (file1)
	file2_lst = bed_to_list (file2)


	logging.info("Gathering information for \"%s\" ..." % file1)
	info1 = bed_info(file1)
	results['A.name'] = info1['Name']
	results['A.interval_count'] = info1['Count']
	results['A.interval_total_size'] = info1['Total_size']
	#results['A.interval_uniq_size'] = info1['Genomic_size']
	results['A.interval_mean_size'] = info1['Mean_size']
	results['A.interval_median_size'] = info1['Median_size']
	results['A.interval_min_size'] = info1['Min_size']
	results['A.interval_max_size'] = info1['Max_size']
	results['A.interval_size_SD'] = info1['STD']
	uniqBase1 = info1['Genomic_size']

	logging.info("Gathering information for \"%s\" ..." % file2)
	info2 = bed_info(file2)
	results['B.name'] = info2['Name']
	results['B.interval_count'] = info2['Count']
	results['B.interval_total_size'] = info2['Total_size']
	#results['B.interval_uniq_size'] = info2['Genomic_size']
	results['B.interval_mean_size'] = info2['Mean_size']
	results['B.interval_median_size'] = info2['Median_size']
	results['B.interval_min_size'] = info2['Min_size']
	results['B.interval_max_size'] = info2['Max_size']
	results['B.interval_size_SD'] = info2['STD']
	uniqBase2 = info2['Genomic_size']

	#calculate overall collocation coef
	logging.debug("Calculating overlapped bases ...")
	overlapBases = bed_overlap_size(file1_lst, file2_lst)

	results['G.size'] = bg_size
	results['A.size'] = uniqBase1
	results['Not_A.size'] = bg_size - uniqBase1
	results['B.size'] = uniqBase2
	results['Not_B.size'] = bg_size - uniqBase2
	results['A_not_B.size'] = uniqBase1 - overlapBases
	results['B_not_A.size'] = uniqBase2 - overlapBases
	results['A_and_B.size'] = overlapBases
	results['A_and_B.exp_size'] = uniqBase1 * uniqBase2/bg_size
	results['A_or_B.size'] = uniqBase1 + uniqBase2 -overlapBases
	results['Neither_A_nor_B.size'] = bg_size - uniqBase1 - uniqBase2 + overlapBases


	results['coef.Collocation'] = ov_coef(uniqBase1, uniqBase2, overlapBases, bg_size)
	results['coef.Jaccard'] = ov_jaccard(uniqBase1, uniqBase2, overlapBases, bg_size)
	results['coef.Dice'] = ov_sd(uniqBase1, uniqBase2, overlapBases, bg_size)
	results['coef.SS'] = ov_ss(uniqBase1, uniqBase2, overlapBases, bg_size)
	results['A_and_B.PMI'] = pmi_value(uniqBase1, uniqBase2, overlapBases, bg_size)
	results['A_and_B.NPMI'] = npmi_value(uniqBase1, uniqBase2, overlapBases, bg_size)


	return pd.Series(data=results)





if __name__ == '__main__':
	pd.set_option('display.float_format', lambda x: '%.4f' % x)
	a = ov_stats(sys.argv[1], sys.argv[2])
	print (a)

