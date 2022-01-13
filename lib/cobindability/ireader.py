
#!/usr/bin/env python
# encoding: utf-8

"""
Read regular, compressed (.gz .bz), remote (http://, https://, ftp://) BED file or BigBed file.
"""


import sys
import bz2
import gzip
from urllib.request import urlopen
from subprocess import Popen, PIPE
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



def bbopen(fname):
	"""
	Open bigBed file. Local and remote bigBed read access is supported.
	"""
	bb = pyBigWig.open(fname)
	chrom_dict = bb.chroms()
	for chr in chrom_dict:
		for start, end, score in bb.entries(chr, 0, chrom_dict[chr]):
			yield(chr + '\t' + str(start) + '\t' + str(end) + '\t' + score)


def nopen(f, mode="rb"):
	"""
	Open regular or compressed BED file.
	"""
	if not isinstance(f, str):
		return f
	if f.startswith("|"):
		p = Popen(f[1:], stdout=PIPE, stdin=PIPE, shell=True)
		if mode[0] == "r": return p.stdout
		return p
	return {"r": sys.stdin, "w": sys.stdout}[mode[0]] if f == "-" \
		else gzip.open(f, mode) if f.endswith((".gz", ".Z", ".z")) \
		else bz2.BZ2File(f, mode) if f.endswith((".bz", ".bz2", ".bzip2")) \
		else urlopen(f) if f.startswith(("http://", "https://","ftp://")) \
		else open(f, mode)


def reader(fname):
	if fname.endswith(('.bb','.bigbed','.bigBed','.BigBed', '.BB',' BIGBED')):
		for l in bbopen(fname):
			yield l
	else:
		for l in nopen(fname):
			yield l.decode('utf8').strip().replace("\r", "")
