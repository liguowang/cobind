#!/usr/bin/env python
# encoding: utf-8

"""
read bigBed files.

"""

import pyBigWig
from cobindability import version

__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "MIT"
__version__ = version.version
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Development"


def bbreader(fname):
    """
    read bigBed file
    """
    bb = pyBigWig.open(fname)
    chrom_dict = bb.chroms()
    for chr in chrom_dict:
        for start, end, score in bb.entries(chr, 0, chrom_dict[chr]):
            yield(chr + '\t' + str(start) + '\t' + str(end) + '\t' + score)
