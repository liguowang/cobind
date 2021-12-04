#!/usr/bin/env python


import sys
#import math
#from bx.bitset import *
#from bx.bitset_builders import *
#import scipy.stats as stats

from cobindability.BED import bedinfo
from cobindability.ovcoef import cal_overlap_coef, cal_pmi
import logging
import argparse

__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="0.0.1"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Development"





def main():
	general_help = "**cobind: cooccurrence, cooccupany, and covariance analysis of genomic regions**"
	bed_help = "Genomic regions in BED (Browser Extensible Data, https://genome.ucsc.edu/FAQ/FAQformat.html#format1), BED-like or BigBed format. The BED-like format includes:\
		'bed3','bed4','bed6','bed12','bedgraph','narrowpeak', 'broadpeak','gappedpeak'. BED and BED-like format can be plain text, compressed (.gz, .z, .bz, .bz2, .bzip2) \
		or remote (http://, https://, ftp://) files. Do not compress BigBed foramt. BigBed file can also be a remote file."
	# sub commands and help.
	commands = {
	'overlap' : "Calculate the overlapping coefficient between two sets of genomic regions.",
	'cooccur' : "Calculate the cooccurrence (mutual information) between two sets of genomic regions.",
	'covary' : "Calculate the covariance (Pearson and Spearman's coefficient) between two sets of genomic regions.",
	'bedinfo' : "Report basic statistics of genomic regions.",
	}

	#create parse
	parser = argparse.ArgumentParser(description=general_help, epilog='', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument('-v', '--version', action='version', version='%s %s' %  ('cobind', __version__))


	# create sub-parser
	sub_parsers = parser.add_subparsers(help='Sub-command description:')
	parser_overlap = sub_parsers.add_parser('overlap', help=commands['overlap'])
	parser_cooccur = sub_parsers.add_parser('cooccur', help=commands['cooccur'])
	parser_covary = sub_parsers.add_parser('covary', help=commands['covary'])
	parser_bedinfo = sub_parsers.add_parser('bedinfo', help=commands['bedinfo'])


	# create the parser for the "overlap" sub-command
	parser_overlap.add_argument("region1", type=str, metavar ="input1.bed",help=bed_help)
	parser_overlap.add_argument("region2", type=str, metavar ="input2.bed",help=bed_help)
	parser_overlap.add_argument('-n', '--ndraws', type=int, dest="iter", default = 50, help="Times of resampling to estimate confidence intervals. Set to '0' to turn off resampling.(default: %(default)d)")
	parser_overlap.add_argument('-s', '--size', type=int, dest="subsize", default = 0.85, help="Size of the subset during resampling. If the original BED file contains 10000 regions, '--size = 0.85' means 8500 regions will be resampled. (default: %(default).2f)")
	parser_overlap.add_argument('-b', '--background', type=int, dest="bgsize", default = 1.4e9, help="The size of the cis-regulatory genomic regions. This is about 1.4Gb For the human genome. (default: %(default)d)")
	parser_overlap.add_argument("-d", "--debug",action="store_true", help="Print detailed information for debugging.")


	# create the parser for the "cooccur" sub-command
	parser_cooccur.add_argument("region1", type=str, metavar ="input1.bed",help=bed_help)
	parser_cooccur.add_argument("region2", type=str, metavar ="input2.bed",help=bed_help)
	parser_cooccur.add_argument('-b', '--background', type=int, dest="bgsize", default = 1.4e9, help="The size of the cis-regulatory genomic regions. This is about 1.4Gb For the human genome. (default: %(default)d)")
	parser_cooccur.add_argument("-d", "--debug",action="store_true", help="Print detailed information for debugging.")

	# create the parser for the "bedinfo" sub-command
	parser_bedinfo.add_argument("bedfile", type=str, metavar ="input.bed",help=bed_help)


	args = parser.parse_args()


	if len(sys.argv)==1:
		parser.print_help(sys.stderr)
		sys.exit(0)
	elif len(sys.argv) >= 2:
		command = sys.argv[1].lower()
		if command =='bedinfo':
			info = bedinfo(args.bedfile)
			print (info)
		elif command == 'overlap':
			if args.debug:
				logging.basicConfig(format = "%(asctime)s [%(levelname)s]  %(message)s",datefmt='%Y-%m-%d %I:%M:%S', level=logging.DEBUG)
			else:
				logging.basicConfig(format = "%(asctime)s [%(levelname)s]  %(message)s",datefmt='%Y-%m-%d %I:%M:%S', level=logging.INFO)
			result = cal_overlap_coef(args.region1, args.region2, n_draws = args.iter, size = args.subsize, bg_size = args.bgsize)
			print (result)
		elif command == 'cooccur':
			if args.debug:
				logging.basicConfig(format = "%(asctime)s [%(levelname)s]  %(message)s",datefmt='%Y-%m-%d %I:%M:%S', level=logging.DEBUG)
			else:
				logging.basicConfig(format = "%(asctime)s [%(levelname)s]  %(message)s",datefmt='%Y-%m-%d %I:%M:%S', level=logging.INFO)
			result = cal_pmi(args.region1, args.region2, bg_size = args.bgsize)
			print (result)

if __name__ == '__main__':
	main()