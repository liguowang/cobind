#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 22:16:41 2022

@author: m102324
"""
import math
import logging
import numpy as np
from cobindability import version


__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "MIT"
__version__ = version.version
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Development"

def ov_coef(x, y, xy):
	"""
	Calculate the overlap coefficient.

	Parameters
	----------
	x : int
		The cardinality of A.
	y : int
		The cardinality of B.
	xy : int
		The cardinality of A AND B.
	Returns
	-------
	float
		Overlap coefficient.
	"""
	if xy > x or xy > y:
		logging.error("Invalid parameters.")
	if min(x,y,xy) < 0:
		logging.error("Invalid parameters.")
	if x == 0 or y == 0 or xy == 0:
		return 0
	else:
		return xy/(math.sqrt(x * y))

def ov_jaccard(x, y, xy):
	"""
	Calculate the Jaccard's coefficient (https://en.wikipedia.org/wiki/Jaccard_index).

	Parameters
	----------
	x : int
		The cardinality of A.
	y : int
		The cardinality of B.
	xy : int
		The cardinality of A AND B.
	Returns
	-------
	float
		Overlap coefficient.
	"""
	if xy > x or xy > y:
		logging.error("Invalid parameters.")
	if min(x,y,xy) < 0:
		logging.error("Invalid parameters.")
	if x == 0 or y == 0:
		return 0
	else:
		return xy/(x + y - xy)

def ov_ss(x, y, xy):
	"""
	Calculate the Szymkiewicz–Simpson coefficient, (https://en.wikipedia.org/wiki/Overlap_coefficient).

	Parameters
	----------
	x : int
		The cardinality of A.
	y : int
		The cardinality of B.
	xy : int
		The cardinality of A AND B.
	Returns
	-------
	float
		Szymkiewicz–Simpson coefficient.
	"""
	if xy > x or xy > y:
		logging.error("Invalid parameters.")
	if min(x,y,xy) < 0:
		logging.error("Invalid parameters.")
	if x == 0 or y == 0:
		return 0
	else:
		return xy/min(x,y)

def ov_sd(x, y, xy):
	"""
	Calculate the Sørensen–Dice coefficient (https://en.wikipedia.org/wiki/Sørensen–Dice_coefficient).
	Sørensen–Dice coefficient also called "Sørensen–Dice index", "Sørensen index" or "Dice's coefficient".

	Parameters
	----------
	x : int
		The cardinality of A.
	y : int
		The cardinality of B.
	xy : int
		The cardinality of A AND B.
	Returns
	-------
	float
		Overlap coefficient.
	"""
	if xy > x or xy > y:
		logging.error("Invalid parameters.")
	if min(x,y,xy) < 0:
		logging.error("Invalid parameters.")
	if x == 0 or y == 0:
		return 0
	else:
		return 2*xy/(x + y)


def pmi_value(x, y, xy, g):
	"""
	Calculate the pointwise mutual information (PMI).

	Parameters
	----------
	x : int
		The cardinality of A.
	y : int
		The cardinality of B.
	xy : int
		The cardinality of A AND B.
	g : int
		The cardinality of background.

	Returns
	-------
	float
		Pointwise mututal information.
	"""
	if g <= 0:
		logging.error("The cardinality of background must be a postive integer.")
	if xy > x or xy > y:
		logging.error("Invalid parameters.")
	if x > g or y > g:
		logging.error("Invalid parameters.")
	if min(x,y,xy) < 0:
		logging.error("Invalid parameters.")
	px = x/g
	py = y/g
	pxy = xy/g

	if pxy == 0:
		return -np.inf
	else:
		return np.log2(pxy) - np.log2(px) - np.log2(py)

def npmi_value(x, y, xy, g):
	"""
	Calculate the normalized pointwise mutual information (NPMI).

	Parameters
	----------
	x : int
		The cardinality of A.
	y : int
		The cardinality of B.
	xy : int
		The cardinality of A AND B.
	g : int
		The cardinality of background.

	Returns
	-------
	float
		Normalized pointwise mutual information.
	"""

	if g <= 0:
		logging.error("The cardinality of background must be a postive integer.")
	if xy > x or xy > y:
		logging.error("Invalid parameters.")
	if x > g or y > g:
		logging.error("Invalid parameters.")
	if min(x,y,xy) < 0:
		logging.error("Invalid parameters.")
	px = x/g
	py = y/g
	pxy = xy/g

	if pxy == 0:
		return -1
	else:
		return np.log2(px*py)/np.log2(pxy) - 1