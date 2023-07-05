#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 12 18:50:00 2021

@author: m102324
"""
import os
import logging
import pyBigWig
import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr, kendalltau
from cobindability import version

__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "MIT"
__version__ = version.version
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Development"


def bigwig_corr(bed, bw1, bw2, outfile, na_label='nan', score_type='mean',
                exact_scores=True, keep_NA=False, top_x=1.0, min_sig=0):
    """
    Calculate Pearson's and Spearman's correlations between two bigWig files

    Parameters
    ----------
    bed : list
        List of genomic regions.
    bw1 : str
        Name of one bigWig file.
    bw2 : str
        Name of another bigWig file.
    outfile : str
        Name of the output TSV file.
    na_label : str, optional
        String representation of missing values. Default: 'nan'
    score_type : str, optional
        Summary statistic score type ('min','mean' or 'max') of a genomic
        region. Default: 'mean'
    exact_scores : bool, optional
        If set, calculate the "exact" summary statistic score rather than
        "zoom-level" score for each genomic region. Default: True
    keep_NA : bool, optional
        If set, all genomic regions will be kept even they do not have summary
        statistic scores in either of the two bigWig files. Default: False
    top_x : float, optional
        Percentage ( if top_x in (0,1]) or number (if top_x > 1) of genomic
        regions used to calculate Pearson, Spearman, Kendall's correlations.
        default: 1.0 (i.e., use all genomic regions to calculate correlation
        coefficients).
    min_sig : float, optional
        Genomic regions with summary statistic equal or less than this
        value will be filtered out. default: 0 (i.e., use all genomic regions
        to calculate correlation coefficients).


    Returns
    -------
    None.

    """
    bw_1 = pyBigWig.open(bw1)
    logging.debug(str(bw_1.header()))
    bw_2 = pyBigWig.open(bw2)
    logging.debug(str(bw_2.header()))

    # check bigwig files
    if bw_1.isBigWig():
        pass
    else:
        logging.error("Not a bigWig file: %s" % bw1)
    if bw_2.isBigWig():
        pass
    else:
        logging.error("Not a bigWig file: %s" % bw2)

    # check chrom IDs in bigwig files
    # all_chroms1 = bw_1.chroms().keys()
    # all_chroms2 = bw_2.chroms().keys()

    names = []
    scores_1 = []
    scores_2 = []
    for (chrom, start, end) in bed:
        region_id = chrom + ':' + str(start) + '-' + str(end)
        score_1 = bw_1.stats(
            chrom, start, end, type=score_type, exact=exact_scores).pop()
        score_2 = bw_2.stats(
            chrom, start, end, type=score_type, exact=exact_scores).pop()
        if (isinstance(score_1, (int, float)) is False):
            if keep_NA:
                score_1 = np.nan
            else:
                continue
        if (isinstance(score_2, (int, float)) is False):
            if keep_NA:
                score_2 = np.nan
            else:
                continue
        names.append(region_id)
        scores_1.append(score_1)
        scores_2.append(score_2)

    bw1_name = os.path.basename(bw1) + '.' + score_type
    bw2_name = os.path.basename(bw2) + '.' + score_type
    df = pd.DataFrame(data={bw1_name: scores_1,
                            bw2_name: scores_2},
                      index=names,
                      dtype=float)

    logging.info("Sort dataframe by summary statistical scores ...")
    df.sort_values(by=[bw1_name, bw1_name],
                   ascending=False,
                   ignore_index=False,
                   inplace=True)

    logging.info("Save dataframe to: \"%s\"" % outfile)
    df.to_csv(outfile,
              na_rep=na_label,
              sep="\t",
              index=True,
              header=True,
              index_label="region_id")

    # Still remove NAs, in order to calculate correlations
    # if keep_NA:
    df = df[df > min_sig]
    df.dropna(axis=0, inplace=True)

    if top_x > 0 and top_x <= 1:
        top_n = int(len(df) * top_x)
        logging.info("Select %d regions ..." % top_n)
        df = df.head(top_n)
    elif top_x > 1:
        top_n = top_x
        logging.info("Select %d regions ..." % top_n)
        df = df.head(top_n)

    (pearson_cor, pearson_p) = pearsonr(
        np.log2(df[bw1_name]), np.log2(df[bw2_name]))
    (spearman_rho, spearman_p) = spearmanr(
        np.log2(df[bw1_name]), np.log2(df[bw2_name]))
    (kendall_tau, kendall_p) = kendalltau(
        np.log2(df[bw1_name]), np.log2(df[bw2_name]))

    return (pd.DataFrame(data={'Pearson_cor:': [pearson_cor, pearson_p],
                               'Spearman_rho:': [spearman_rho, spearman_p],
                               'Kendall_tau:': [kendall_tau, kendall_p]},
                         index=['Correlation', 'P-value']))
