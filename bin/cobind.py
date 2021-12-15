#!/usr/bin/env python

import sys

from cobindability.BED import bedinfo,overlap_bed,bedtolist
from cobindability.ovcoef import cal_overlap_coef, cal_pmi
from cobindability.bw import bigwig_corr
import logging
import argparse

__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="0.0.2"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Development"

def main():
	general_help = "**cobind: collocation analyses of genomic regions**"
	bed_help = "Genomic regions in BED, BED-like or bigBed format. The BED-like format includes:\
		'bed3', 'bed4', 'bed6', 'bed12', 'bedgraph', 'narrowpeak', 'broadpeak', 'gappedpeak'. BED and BED-like format can be plain text, compressed (.gz, .z, .bz, .bz2, .bzip2) \
		or remote (http://, https://, ftp://) files. Do not compress BigBed foramt. BigBed file can also be a remote file."
	# sub commands and help.
	commands = {
	'overlap' : "Calculate the overlapping coefficient between two sets of genomic regions.",
	'pmi' : "Calculate the pointwise mutual information (PMI) between two sets of genomic regions.",
	'cooccur' : "Calculate the co-occurrence between two sets of genomic regions.",
	'covary' : "Calculate the covariance (Pearson, Spearman and Kendall coefficients) of binding intensities between two sets of genomic regions.",
	'bedinfo' : "Report basic statistics of genomic regions.",
	}

	#create parse
	parser = argparse.ArgumentParser(description=general_help, epilog='', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument('-v', '--version', action='version', version='%s %s' %  ('cobind', __version__))


	# create sub-parser
	sub_parsers = parser.add_subparsers(help='Sub-command description:')
	parser_overlap = sub_parsers.add_parser('overlap', help=commands['overlap'])
	parser_pmi = sub_parsers.add_parser('pmi', help=commands['pmi'])
	parser_cooccur = sub_parsers.add_parser('cooccur', help=commands['cooccur'])
	parser_covary = sub_parsers.add_parser('covary', help=commands['covary'])
	parser_bedinfo = sub_parsers.add_parser('bedinfo', help=commands['bedinfo'])

	# create the parser for the "overlap" sub-command
	parser_overlap.add_argument("bed1", type=str, metavar ="input_A.bed",help=bed_help)
	parser_overlap.add_argument("bed2", type=str, metavar ="input_B.bed",help=bed_help)
	parser_overlap.add_argument('-n', '--ndraws', type=int, dest="iter", default = 50, help="Times of resampling to estimate confidence intervals. Set to '0' to turn off resampling.(default: %(default)d)")
	parser_overlap.add_argument('-s', '--size', type=int, dest="subsize", default = 0.85, help="Size of the subset during resampling. If the original BED file contains 10000 regions, '--size = 0.85' means 8500 regions will be resampled. (default: %(default).2f)")
	parser_overlap.add_argument('-b', '--background', type=int, dest="bgsize", default = 1.4e9, help="The size of the cis-regulatory genomic regions. This is about 1.4Gb For the human genome. (default: %(default)d)")
	parser_overlap.add_argument("-d", "--debug",action="store_true", help="Print detailed information for debugging.")

	# create the parser for the "pmi" sub-command
	parser_pmi.add_argument("bed1", type=str, metavar ="input_A.bed",help=bed_help)
	parser_pmi.add_argument("bed2", type=str, metavar ="input_B.bed",help=bed_help)
	parser_pmi.add_argument('-b', '--background', type=int, dest="bgsize", default = 1.4e9, help="The size of the cis-regulatory genomic regions. This is about 1.4Gb For the human genome. (default: %(default)d)")
	parser_pmi.add_argument("-d", "--debug",action="store_true", help="Print detailed information for debugging.")

	# create the parser for the "cooccur" sub-command
	parser_cooccur.add_argument("bed1", type=str, metavar ="input_A.bed",help=bed_help)
	parser_cooccur.add_argument("bed2", type=str, metavar ="input_B.bed",help=bed_help)
	parser_cooccur.add_argument('-b', '--background', type=int, dest="bgsize", default = 1.4e9, help="The size of the cis-regulatory genomic regions. This is about 1.4Gb For the human genome. (default: %(default)d)")
	parser_cooccur.add_argument("-d", "--debug",action="store_true", help="Print detailed information for debugging.")

	# create the parser for the "covary" sub-command
	parser_covary.add_argument("bed1", type=str, metavar="input_A.bed",help=bed_help)
	parser_covary.add_argument("bw1", type=str, metavar="input_A.bw",help="Input bigWig file matched to 'input_A.bed'. BigWig file can be local or remote. Note: the chromosome IDs must be consistent between BED and bigWig files.")
	parser_covary.add_argument("bed2", type=str, metavar="input_B.bed",help=bed_help)
	parser_covary.add_argument("bw2", type=str, metavar="input_B.bw",help="Input bigWig file matched to 'input_B.bed'. BigWig file can be local or remote.  Note: the chromosome IDs must be consistent between BED and bigWig files.")
	parser_covary.add_argument("output", type=str,metavar="output_prefix", help="Prefix of output files. Three files will be generated: \"output_prefix_bedA_unique.tsv\" (input_A.bed specific regions and their bigWig scores), \"output_prefix_bedB_unique.tsv\" (input_B.bed specific regions and their bigWig scores), and \"output_prefix_common.tsv\"(input_A.bed and input_B.bed overlapped regions and their bigWig scores).")
	parser_covary.add_argument("--na", type=str, dest="na_label", default = 'nan', help="Symbols used to represent the missing values. (default: %(default)s)")
	parser_covary.add_argument('--type', type=str, dest="score_type", choices=['mean', 'min', 'max'],  default = 'mean', help="Summary statistic score type ('min','mean' or 'max') of a genomic region. (default: %(default)s)")
	parser_covary.add_argument('--topx', type=float, dest="top_X",  default = 0.1, help="Fraction/percentage (if top_X in (0,1]) or number (if top_X > 1) of genomic regions used to calculate Pearson, Spearman, Kendall's correlations. If TOP_X == 1, all (i.e., 100%%) genomic regions will be used to calculate correlations. (default: %(default)s)")
	parser_covary.add_argument("--exact", dest="exact", action="store_true", help="If set, calculate the \"exact\" summary statistic score rather than \"zoom-level\" score for each genomic region.")
	parser_covary.add_argument("--nosort", dest="nosort", action="store_true", help="If set, will NOT sort the summary statistical scores.")
	parser_covary.add_argument("--keepna", dest="keepna", action="store_true", help="If set, a genomic region will be kept even it does not have summary statistical score in either of the two bigWig files. This flag only affects the output .tsv files.")
	parser_covary.add_argument("-d", "--debug",action="store_true", help="Print detailed information for debugging.")

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
			result = cal_overlap_coef(args.bed1, args.bed2, n_draws = args.iter, size = args.subsize, bg_size = args.bgsize)
			print (result)
		elif command == 'pmi':
			if args.debug:
				logging.basicConfig(format = "%(asctime)s [%(levelname)s]  %(message)s",datefmt='%Y-%m-%d %I:%M:%S', level=logging.DEBUG)
			else:
				logging.basicConfig(format = "%(asctime)s [%(levelname)s]  %(message)s",datefmt='%Y-%m-%d %I:%M:%S', level=logging.INFO)
			result = cal_pmi(args.bed1, args.bed2, bg_size = args.bgsize)
			print (result)
		elif command == 'covary':
			if args.debug:
				logging.basicConfig(format = "%(asctime)s [%(levelname)s]  %(message)s",datefmt='%Y-%m-%d %I:%M:%S', level=logging.DEBUG)
			else:
				logging.basicConfig(format = "%(asctime)s [%(levelname)s]  %(message)s",datefmt='%Y-%m-%d %I:%M:%S', level=logging.INFO)

			a,b,c = overlap_bed(args.bed1, args.bed2)
			logging.info("Calculate summay statistics for overlapped regions ...")
			c_corr = bigwig_corr(bed = c, bw1 = args.bw1, bw2 = args.bw2, outfile = args.output + '_common.tsv', na_label = args.na_label, score_type = args.score_type, exact_scores = args.exact, keep_NA = args.keepna, no_sort = args.nosort, top_x = args.top_X)
			print (c_corr.T)

			logging.info("Calculate summay statistics for \"%s\" unique regions ..." % args.bed1)
			a_corr = bigwig_corr(bed = a, bw1 = args.bw1, bw2 = args.bw2, outfile = args.output + '_bedA_unique.tsv', na_label = args.na_label, score_type = args.score_type, exact_scores = args.exact,  keep_NA = args.keepna, no_sort = args.nosort, top_x = args.top_X)
			print (a_corr.T)

			logging.info("Calculate summay statistics for \"%s\" unique regions ..." % args.bed2)
			b_corr = bigwig_corr(bed = b, bw1 = args.bw1, bw2 = args.bw2, outfile = args.output + '_bedB_unique.tsv', na_label = args.na_label, score_type = args.score_type, exact_scores = args.exact,  keep_NA = args.keepna, no_sort = args.nosort, top_x = args.top_X)
			print (b_corr.T)
		elif command == 'cooccur':
			pass
if __name__ == '__main__':
	main()
