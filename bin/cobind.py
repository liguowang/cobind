#!/usr/bin/env python

import sys
import logging
import argparse
import pandas as pd
from cobindability.BED import compare_bed, peakwise_ovcoef
from cobindability.BED import cooccur_peak, srog_peak
from cobindability.bw import bigwig_corr
from cobindability.ovstat import ov_stats
from cobindability import version
from cobindability.ovbootstrap import bootstrap_coef, bootstrap_npmi
from cobindability.coefcal import ov_coef, ov_jaccard, ov_ss, ov_sd
from cobindability.coefcal import pmi_value, npmi_value
from cobindability.utils import config_log, cal_zscores

__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "MIT"
__version__ = version.version
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Production"


def main():
    pd.set_option('display.float_format', lambda x: '%.4f' % x)
    general_help = "**cobind: collocation analyses of genomic regions**"
    bed_help = "Genomic regions in BED, BED-like or bigBed format. The \
        BED-like format includes:'bed3', 'bed4', 'bed6', 'bed12', 'bedgraph', \
        'narrowpeak', 'broadpeak', 'gappedpeak'. BED and BED-like format can \
        be plain text, compressed (.gz, .z, .bz, .bz2, .bzip2) or remote \
        (http://, https://, ftp://) files. Do not compress BigBed foramt. \
        BigBed file can also be a remote file."
    log_help = "\
        This file is used to save the log information. By default, if no file \
        is specified (None), the log information will be printed to the \
        screen."
    nameA_help = "\
        Name to represent 1st set of genomic interval. If not \
        specified (None), the file name (\"input_A.bed\") will be used."
    nameB_help = "\
        Name to represent the 2nd set of genomic interval. If not \
        specified (None), the file name (\"input_B.bed\") will be used."
    # sub commands and help.
    commands = {
        'overlap': "Calculate the collocation coefficient (C) between two \
            sets of genomic regions. C = |A and B| / (|A|*|B|)**0.5",
        'jaccard': "Calculate the Jaccard similarity coefficient (J) between \
            two sets of genomic regions. J = |A and B| / |A or B|",
        'dice': "Calculate the Sørensen–Dice coefficient (SD) between two \
            sets of genomic regions. SD = 2*|A and B| / (|A| + |B|)",
        'simpson': "Calculate the Szymkiewicz–Simpson coefficient (SS) \
            between two sets of genomic regions. \
            SS = |A and B| / min(|A|, |B|)",
        'pmi': "Calculate the pointwise mutual information (PMI) between two \
            sets of genomic regions. PMI = log(p(|A and B|)) - log(p(|A|)) - \
            log(p(|B|))",
        'npmi': "Calculate the normalized pointwise mutual information (NPMI) \
            between two sets of genomic regions. NPMI = log(p(|A|)*p(|B|)) / \
            log(p(|A and B|)) - 1",
        'cooccur': "Evaluate if two sets of genomic regions are significantly \
            co-occurred in given background regions.",
        'covary': "Calculate the covariance (Pearson, Spearman and Kendall \
            coefficients) of binding intensities between two sets of genomic \
            regions.",
        'srog': "Report the code of Spatial Relation Of Genomic (SROG) \
            regions. SROG codes include 'disjoint','touch','equal','overlap',\
            'contain', 'within'.",
        'stat': "Wrapper function. Report basic statistics of genomic \
            regions, and calculate overlapping measurements (including \"C\", \
            \"J\", \"SD\", \"SS\", \"PMI\", \"NPMI\"), without bootstrap \
            resampling or generating peakwise measurements.",
        'zscore': "Calculate Z-score of six overlapping measurements \
            (inlcuding \"C\", \"J\", \"SD\", \"SS\", \"PMI\", \"NPMI\"), to \
            provide an overall measurement of the collocation strength."
    }

    # create parse
    parser = argparse.ArgumentParser(
        description=general_help, epilog='',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%s %s' % ('cobind', __version__)
        )

    # create sub-parser
    sub_parsers = parser.add_subparsers(
        help='Sub-command description:')
    parser_overlap = sub_parsers.add_parser(
        'overlap', help=commands['overlap'])
    parser_jaccard = sub_parsers.add_parser(
        'jaccard', help=commands['jaccard'])
    parser_dice = sub_parsers.add_parser(
        'dice', help=commands['dice'])
    parser_simpson = sub_parsers.add_parser(
        'simpson', help=commands['simpson'])
    parser_pmi = sub_parsers.add_parser(
        'pmi', help=commands['pmi'])
    parser_npmi = sub_parsers.add_parser(
        'npmi', help=commands['npmi'])
    parser_cooccur = sub_parsers.add_parser(
        'cooccur', help=commands['cooccur'])
    parser_covary = sub_parsers.add_parser(
        'covary', help=commands['covary'])
    parser_srog = sub_parsers.add_parser(
        'srog', help=commands['srog'])
    parser_stat = sub_parsers.add_parser(
        'stat', help=commands['stat'])
    parser_zscore = sub_parsers.add_parser(
        'zscore', help=commands['zscore'])

    # create the parser for the "overlap" sub-command
    parser_overlap.add_argument(
        "bed1", type=str, metavar="input_A.bed", help=bed_help)
    parser_overlap.add_argument(
        "bed2", type=str, metavar="input_B.bed", help=bed_help)
    parser_overlap.add_argument(
        '--nameA', type=str, default=None, help=nameA_help)
    parser_overlap.add_argument(
        '--nameB', type=str, default=None, help=nameB_help)
    parser_overlap.add_argument(
        '-n', '--ndraws', type=int, dest="iter", default=20,
        help="Times of resampling to estimate confidence intervals. \
            Set to '0' to turn off resampling. For the resampling process \
            to work properly, overlapped intervals in each bed file must be \
            merged. (default: %(default)d)")
    parser_overlap.add_argument(
        '-f', '--fraction', type=int, dest="subsample", default=0.75,
        help="Resampling fraction. (default: %(default).2f)")
    parser_overlap.add_argument(
        '-b', '--background', type=int, dest="bgsize", default=1.4e9,
        help="The size of the cis-regulatory genomic regions. This is about \
            1.4Gb For the human genome. (default: %(default)d)")
    parser_overlap.add_argument(
        "-o", "--save", action="store_true",
        help="If set, will save peak-wise coefficients to files \
            (\"input_A_peakwise_scores.tsv\" and \
             \"input_B_peakwise_scores.tsv\").")
    parser_overlap.add_argument(
        "-l", "--log", type=str, metavar="log_file", default=None,
        help=log_help)
    parser_overlap.add_argument(
        "-d", "--debug", action="store_true",
        help="Print detailed information for debugging.")

    # create the parser for the "jaccard" sub-command
    parser_jaccard.add_argument(
        "bed1", type=str, metavar="input_A.bed", help=bed_help)
    parser_jaccard.add_argument(
        "bed2", type=str, metavar="input_B.bed", help=bed_help)
    parser_jaccard.add_argument(
        '--nameA', type=str, default=None, help=nameA_help)
    parser_jaccard.add_argument(
        '--nameB', type=str, default=None, help=nameB_help)
    parser_jaccard.add_argument(
        '-n', '--ndraws', type=int, dest="iter", default=20,
        help="Times of resampling to estimate confidence intervals. Set to \
            '0' to turn off resampling. For the resampling process to work \
            properly, overlapped intervals in each bed file must be merged. \
            (default: %(default)d)")
    parser_jaccard.add_argument(
        '-f', '--fraction', type=int, dest="subsample", default=0.75,
        help="Resampling fraction. (default: %(default).2f)")
    parser_jaccard.add_argument(
        '-b', '--background', type=int, dest="bgsize", default=1.4e9,
        help="The size of the cis-regulatory genomic regions. This is about \
            1.4Gb For the human genome. (default: %(default)d)")
    parser_jaccard.add_argument(
        "-o", "--save", action="store_true",
        help="If set, will save peak-wise coefficients to files \
            (\"input_A_peakwise_scores.tsv\" and \
             \"input_B_peakwise_scores.tsv\").")
    parser_jaccard.add_argument(
        "-l", "--log", type=str, metavar="log_file", default=None,
        help=log_help)
    parser_jaccard.add_argument(
        "-d", "--debug", action="store_true",
        help="Print detailed information for debugging.")

    # create the parser for the "dice" sub-command
    parser_dice.add_argument(
        "bed1", type=str, metavar="input_A.bed", help=bed_help)
    parser_dice.add_argument(
        "bed2", type=str, metavar="input_B.bed", help=bed_help)
    parser_dice.add_argument(
        '--nameA', type=str, default=None, help=nameA_help)
    parser_dice.add_argument(
        '--nameB', type=str, default=None, help=nameB_help)
    parser_dice.add_argument(
        '-n', '--ndraws', type=int, dest="iter", default=20,
        help="Times of resampling to estimate confidence intervals. Set to \
            '0' to turn off resampling. For the resampling process to work \
            properly, overlapped intervals in each bed file must be  merged. \
            (default: %(default)d)")
    parser_dice.add_argument(
        '-f', '--fraction', type=int, dest="subsample", default=0.75,
        help="Resampling fraction. (default: %(default).2f)")
    parser_dice.add_argument(
        '-b', '--background', type=int, dest="bgsize", default=1.4e9,
        help="The size of the cis-regulatory genomic regions. This is about \
            1.4Gb For the human genome. (default: %(default)d)")
    parser_dice.add_argument(
        "-o", "--save", action="store_true",
        help="If set, will save peak-wise coefficients to files \
            (\"input_A_peakwise_scores.tsv\" and \
             \"input_B_peakwise_scores.tsv\").")
    parser_dice.add_argument(
        "-l", "--log", type=str, metavar="log_file", default=None,
        help=log_help)
    parser_dice.add_argument(
        "-d", "--debug", action="store_true",
        help="Print detailed information for debugging.")

    # create the parser for the "simpson" sub-command
    parser_simpson.add_argument(
        "bed1", type=str, metavar="input_A.bed", help=bed_help)
    parser_simpson.add_argument(
        "bed2", type=str, metavar="input_B.bed", help=bed_help)
    parser_simpson.add_argument(
        '--nameA', type=str, default=None, help=nameA_help)
    parser_simpson.add_argument(
        '--nameB', type=str, default=None, help=nameB_help)
    parser_simpson.add_argument(
        '-n', '--ndraws', type=int, dest="iter", default=20,
        help="Times of resampling to estimate confidence intervals. Set to \
            '0' to turn off resampling. For the resampling process to work \
            properly, overlapped intervals in each bed file must be merged. \
            (default: %(default)d)")
    parser_simpson.add_argument(
        '-f', '--fraction', type=int, dest="subsample", default=0.75,
        help="Resampling fraction. (default: %(default).2f)")
    parser_simpson.add_argument(
        '-b', '--background', type=int, dest="bgsize", default=1.4e9,
        help="The size of the cis-regulatory genomic regions. This is about \
            1.4Gb For the human genome. (default: %(default)d)")
    parser_simpson.add_argument(
        "-o", "--save", action="store_true",
        help="If set, will save peak-wise coefficients to files \
            (\"input_A_peakwise_scores.tsv\" and \
             \"input_B_peakwise_scores.tsv\").")
    parser_simpson.add_argument(
        "-l", "--log", type=str, metavar="log_file", default=None,
        help=log_help)
    parser_simpson.add_argument(
        "-d", "--debug", action="store_true",
        help="Print detailed information for debugging.")

    # create the parser for the "pmi" sub-command
    parser_pmi.add_argument(
        "bed1", type=str, metavar="input_A.bed", help=bed_help)
    parser_pmi.add_argument(
        "bed2", type=str, metavar="input_B.bed", help=bed_help)
    parser_pmi.add_argument(
        '--nameA', type=str, default=None, help=nameA_help)
    parser_pmi.add_argument(
        '--nameB', type=str, default=None, help=nameB_help)
    parser_pmi.add_argument(
        '-n', '--ndraws', type=int, dest="iter", default=20,
        help="Times of resampling to estimate confidence intervals. Set to \
            '0' to turn off resampling. For the resampling process to work \
            properly, overlapped intervals in each bed file must be merged. \
            (default: %(default)d)")
    parser_pmi.add_argument(
        '-f', '--fraction', type=int, dest="subsample", default=0.75,
        help="Resampling fraction. (default: %(default).2f)")
    parser_pmi.add_argument(
        '-b', '--background', type=int, dest="bgsize", default=1.4e9,
        help="The size of the cis-regulatory genomic regions. This is about \
            1.4Gb For the human genome. (default: %(default)d)")
    parser_pmi.add_argument(
        "-o", "--save", action="store_true",
        help="If set, will save peak-wise coefficients to files \
            (\"input_A_peakwise_scores.tsv\" and \
             \"input_B_peakwise_scores.tsv\").")
    parser_pmi.add_argument(
        "-l", "--log", type=str, metavar="log_file", default=None,
        help=log_help)
    parser_pmi.add_argument(
        "-d", "--debug", action="store_true",
        help="Print detailed information for debugging.")

    # create the parser for the "npmi" sub-command
    parser_npmi.add_argument(
        "bed1", type=str, metavar="input_A.bed", help=bed_help)
    parser_npmi.add_argument(
        "bed2", type=str, metavar="input_B.bed", help=bed_help)
    parser_npmi.add_argument(
        '--nameA', type=str, default=None, help=nameA_help)
    parser_npmi.add_argument(
        '--nameB', type=str, default=None, help=nameB_help)
    parser_npmi.add_argument(
        '-n', '--ndraws', type=int, dest="iter", default=20,
        help="Times of resampling to estimate confidence intervals. Set to \
            '0' to turn off resampling. For the resampling process to work \
            properly, overlapped intervals in each bed file must be merged. \
            (default: %(default)d)")
    parser_npmi.add_argument(
        '-f', '--fraction', type=int, dest="subsample", default=0.75,
        help="Resampling fraction. (default: %(default).2f)")
    parser_npmi.add_argument(
        '-b', '--background', type=int, dest="bgsize", default=1.4e9,
        help="The size of the cis-regulatory genomic regions. This is about \
            1.4Gb For the human genome. (default: %(default)d)")
    parser_npmi.add_argument(
        "-o", "--save", action="store_true",
        help="If set, will save peak-wise coefficients to files \
            (\"input_A_peakwise_scores.tsv\" and \
             \"input_B_peakwise_scores.tsv\").")
    parser_npmi.add_argument(
        "-l", "--log", type=str, metavar="log_file", default=None,
        help=log_help)
    parser_npmi.add_argument(
        "-d", "--debug", action="store_true",
        help="Print detailed information for debugging.")

    # create the parser for the "cooccur" sub-command
    parser_cooccur.add_argument(
        "bed1", type=str, metavar="input_A.bed", help=bed_help)
    parser_cooccur.add_argument(
        "bed2", type=str, metavar="input_B.bed", help=bed_help)
    parser_cooccur.add_argument(
        '--nameA', type=str, default=None,
        help="Name to represent 1st set of genomic interval. If not specified\
            \"A\" will be used.")
    parser_cooccur.add_argument(
        '--nameB', type=str, default=None,
        help="Name to represent 2nd set of genomic interval. If not specified\
            \"B\" will be used.")
    parser_cooccur.add_argument(
        "bed3", type=str, metavar="background.bed",
        help="Genomic regions as the background (e.g., all promoters, \
            all enhancers).")
    parser_cooccur.add_argument(
        "output", type=str, metavar="output.tsv",
        help="For each genomic region in the \"background.bed\" file, add \
            another column indicating if this region is \"input_A specific \
            (i.e., A+B-)\", \"input_B specific (i.e., A-B+)\", \"co-occur \
            (i.e., A+B+)\" or \"neither (i.e, A-B-)\". ")
    parser_cooccur.add_argument(
        '--ncut', type=int, dest="n_cut",  default=1,
        help="The minimum overlap size. (default: %(default)d)")
    parser_cooccur.add_argument(
        '--pcut', type=float, dest="p_cut", default=0.0,
        help="The minimum overlap percentage. (default: %(default)f)")
    parser_cooccur.add_argument(
        "-l", "--log", type=str, metavar="log_file", default=None,
        help=log_help)
    parser_cooccur.add_argument(
        "-d", "--debug", action="store_true",
        help="Print detailed information for debugging.")

    # create the parser for the "covary" sub-command
    parser_covary.add_argument(
        "bed1", type=str, metavar="input_A.bed", help=bed_help)
    parser_covary.add_argument(
        "bw1", type=str, metavar="input_A.bw",
        help="Input bigWig file matched to 'input_A.bed'. BigWig file can be \
            local or remote. Note: the chromosome IDs must be consistent \
            between BED and bigWig files.")
    parser_covary.add_argument(
        "bed2", type=str, metavar="input_B.bed", help=bed_help)
    parser_covary.add_argument(
        "bw2", type=str, metavar="input_B.bw",
        help="Input bigWig file matched to 'input_B.bed'. BigWig file can be \
            local or remote.  Note: the chromosome IDs must be consistent \
            between BED and bigWig files.")
    parser_covary.add_argument(
        "output", type=str, metavar="output_prefix",
        help="Prefix of output files. Three files will be generated: \
        \"output_prefix_bedA_unique.tsv\" (input_A.bed specific regions and \
        their bigWig scores), \"output_prefix_bedB_unique.tsv\" (input_B.bed \
        specific regions and their bigWig scores), and \
        \"output_prefix_common.tsv\"(input_A.bed and input_B.bed overlapped \
        regions and their bigWig scores).")
    parser_covary.add_argument(
        '--nameA', type=str, default=None,
        help="Name of the 1st set of genomic interval, if not proviced, \
            \"bedA\" will be used. Only affects the name of output file.")
    parser_covary.add_argument(
        '--nameB', type=str, default=None,
        help="Name of the 2nd set of genomic interval, if not proviced, \
            \"bedB\" will be used. Only affects the name of output file.")
    parser_covary.add_argument(
        "--na", type=str, dest="na_label", default='nan',
        help="Symbols used to represent the missing values. \
            (default: %(default)s)")
    parser_covary.add_argument(
        '--type', type=str, dest="score_type",
        choices=['mean', 'min', 'max'], default='mean',
        help="Summary statistic score type ('min','mean' or 'max') of a \
            genomic region. (default: %(default)s)")
    parser_covary.add_argument(
        '--topx', type=float, dest="top_X", default=1.0,
        help="Fraction (if 0 < top_X <= 1) or number (if top_X > 1) of \
            genomic regions used to calculate Pearson, Spearman, Kendall's \
            correlations. If TOP_X == 1 (i.e., 100%%), all the genomic \
            regions will be used to calculate correlations. \
            (default: %(default)s)")
    parser_covary.add_argument(
        '--min_sig', type=float, dest="min_signal", default=0,
        help="Genomic region with summary statistic score <= this will be \
            removed. (default: %(default)s)")
    parser_covary.add_argument(
        "--exact", dest="exact", action="store_true",
        help="If set, calculate the \"exact\" summary statistic score rather \
            than \"zoom-level\" score for each genomic region.")
    parser_covary.add_argument(
        "--keepna", dest="keepna", action="store_true",
        help="If set, a genomic region will be kept even it does not have \
            summary statistical score in either of the two bigWig files. \
            This flag only affects the output TSV files.")
    parser_covary.add_argument(
        "-l", "--log", type=str, metavar="log_file", default=None,
        help=log_help)
    parser_covary.add_argument(
        "-d", "--debug", action="store_true",
        help="Print detailed information for debugging.")

    # create the parser for the "srog" sub-command
    parser_srog.add_argument(
        "bed1", type=str, metavar="input_A.bed",
        help="Genomic regions in BED, BED-like or bigBed format. If 'name' \
            (the 4th column) is not provided, the default name is \
            \"chrom:start-end\". If strand (the 6th column) is not provided, \
            the default strand is \"+\".")
    parser_srog.add_argument(
        "bed2", type=str, metavar="input_B.bed",
        help="Genomic regions in BED, BED-like or bigBed format. If 'name' \
            (the 4th column) is not provided, the default name is \
            \"chrom:start-end\". If strand (the 6th column) is not provided, \
            the default strand is \"+\". ")
    parser_srog.add_argument(
        "output", type=str, metavar="output.tsv",
        help="Generate spatial relation code (disjoint, touch, equal, \
            overlap, contain, within) for each genomic interval in \
            \"input_A.bed\".")
    parser_srog.add_argument(
        '--dist', type=int, dest="max_dist", default=250000000,
        help="When intervals are disjoint, find the closest up- and \
            down-stream intervals that are no further than `max_dist` away. \
            default: %(default)d)")
    parser_srog.add_argument(
        "-l", "--log", type=str, metavar="log_file", default=None,
        help=log_help)
    parser_srog.add_argument(
        "-d", "--debug", action="store_true",
        help="Print detailed information for debugging.")

    # create the parser for the "stat" sub-command
    parser_stat.add_argument(
        "bed1", type=str, metavar="input_A.bed", help=bed_help)
    parser_stat.add_argument(
        "bed2", type=str, metavar="input_B.bed", help=bed_help)
    parser_stat.add_argument(
        '--nameA', type=str, default=None, help=nameA_help)
    parser_stat.add_argument(
        '--nameB', type=str, default=None, help=nameB_help)
    parser_stat.add_argument(
        '-b', '--background', type=int, dest="bgsize",
        default=1.4e9, help="The size of the cis-regulatory genomic \
            regions. This is about 1.4Gb For the human genome. \
            (default: %(default)d)")
    parser_stat.add_argument(
        "-l", "--log", type=str, metavar="log_file", default=None,
        help=log_help)
    parser_stat.add_argument(
        "-d", "--debug", action="store_true",
        help="Print detailed information for debugging.")

    # create the parser for the "zscore" sub-command
    parser_zscore.add_argument(
        "input", type=str, metavar="input_file.tsv",
        help="Input dataframe with row names and column names. Must separate \
            different columns with tab. If  \"C\", \"J\", \"SD\", \"SS\", \
            \"PMI\", \"NPMI\" are used as the column names, only these six \
            columns will be used to calculate the Z-score, otherwise, all \
            numerical columns in the dataframe will be used.")
    parser_zscore.add_argument(
        "output", type=str, metavar="output_file.tsv",
        help="Output dataframe with Z-scores as the last column.")
    parser_zscore.add_argument(
        "-l", "--log", type=str, metavar="log_file", default=None,
        help=log_help)

    parser_zscore.add_argument(
        "-d", "--debug", action="store_true",
        help="Print detailed information for debugging.")

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)
    elif len(sys.argv) >= 2:
        command = sys.argv[1]
        if command == 'stat':
            config_log(switch=args.debug, logfile=args.log)
            info = ov_stats(args.bed1, args.bed2,
                            name1=args.nameA,
                            name2=args.nameB,
                            bg_size=args.bgsize)
            print(info)

        elif command == 'overlap':
            config_log(switch=args.debug, logfile=args.log)
            logging.info("Calculate collocation coefficient (overall) ...")
            result = bootstrap_coef(args.bed1,
                                    args.bed2,
                                    name1=args.nameA,
                                    name2=args.nameB,
                                    size_factor=1/args.subsample,
                                    score_func=ov_coef,
                                    n_draws=args.iter,
                                    fraction=args.subsample,
                                    bg_size=args.bgsize)
            print(result)
            if args.save:
                logging.info(
                    "Calculate collocation coefficient (peak-wise) ...")
                peakwise_ovcoef(args.bed1,
                                args.bed2,
                                name1=args.nameA,
                                name2=args.nameB,
                                score_func=ov_coef,
                                g=args.bgsize,
                                na_label='NA')

        elif command == 'jaccard':
            config_log(switch=args.debug, logfile=args.log)
            logging.info("Calculate Jaccard coefficient (overall) ...")
            result = bootstrap_coef(args.bed1,
                                    args.bed2,
                                    name1=args.nameA,
                                    name2=args.nameB,
                                    size_factor=1/args.subsample,
                                    score_func=ov_jaccard,
                                    n_draws=args.iter,
                                    fraction=args.subsample,
                                    bg_size=args.bgsize)
            print(result)
            if args.save:
                logging.info("Calculate Jaccard coefficient (peakwise) ...")
                peakwise_ovcoef(args.bed1,
                                args.bed2,
                                name1=args.nameA,
                                name2=args.nameB,
                                score_func=ov_jaccard,
                                g=args.bgsize,
                                na_label='NA')

        elif command == 'dice':
            config_log(switch=args.debug, logfile=args.log)
            logging.info("Calculate Sørensen–Dice coefficient (overall) ...")
            result = bootstrap_coef(args.bed1,
                                    args.bed2,
                                    name1=args.nameA,
                                    name2=args.nameB,
                                    size_factor=1/args.subsample,
                                    score_func=ov_sd, n_draws=args.iter,
                                    fraction=args.subsample,
                                    bg_size=args.bgsize)
            print(result)
            if args.save:
                logging.info(
                    "Calculate Sørensen–Dice coefficient (peakwise) ...")
                peakwise_ovcoef(args.bed1,
                                args.bed2,
                                name1=args.nameA,
                                name2=args.nameB,
                                score_func=ov_sd,
                                g=args.bgsize,
                                na_label='NA')

        elif command == 'simpson':
            config_log(switch=args.debug, logfile=args.log)
            logging.info(
                "Calculate Szymkiewicz–Simpson coefficient (overall) ...")
            result = bootstrap_coef(args.bed1,
                                    args.bed2,
                                    name1=args.nameA,
                                    name2=args.nameB,
                                    size_factor=1/args.subsample,
                                    score_func=ov_ss,
                                    n_draws=args.iter,
                                    fraction=args.subsample,
                                    bg_size=args.bgsize)
            print(result)
            if args.save:
                logging.info(
                    "Calculate Szymkiewicz–Simpson coefficient (peakwise) ...")
                peakwise_ovcoef(args.bed1,
                                args.bed2,
                                name1=args.nameA,
                                name2=args.nameB,
                                score_func=ov_ss,
                                g=args.bgsize,
                                na_label='NA')

        elif command == 'pmi':
            config_log(switch=args.debug, logfile=args.log)
            logging.info(
                "Calculate the pointwise mutual information (PMI) ...")
            # result = cal_pmi(args.bed1, args.bed2, bg_size=args.bgsize)
            result = bootstrap_coef(args.bed1,
                                    args.bed2,
                                    name1=args.nameA,
                                    name2=args.nameB,
                                    size_factor=1,
                                    score_func=pmi_value,
                                    n_draws=args.iter,
                                    fraction=args.subsample,
                                    bg_size=args.bgsize)
            print(result)
            if args.save:
                peakwise_ovcoef(args.bed1,
                                args.bed2,
                                name1=args.nameA,
                                name2=args.nameB,
                                score_func=pmi_value,
                                g=args.bgsize,
                                na_label='NA')

        elif command == 'npmi':
            config_log(switch=args.debug, logfile=args.log)
            logging.info(
                "Calculate the normalized pointwise mutual information (NPMI)")
            # result=cal_pmi(args.bed1, args.bed2, bg_size=args.bgsize)
            result = bootstrap_npmi(args.bed1,
                                    args.bed2,
                                    name1=args.nameA,
                                    name2=args.nameB,
                                    size_factor=1/args.subsample,
                                    score_func=npmi_value,
                                    n_draws=args.iter,
                                    fraction=args.subsample,
                                    bg_size=args.bgsize)
            print(result)
            if args.save:
                peakwise_ovcoef(args.bed1,
                                args.bed2,
                                name1=args.nameA,
                                name2=args.nameB,
                                score_func=npmi_value,
                                g=args.bgsize,
                                na_label='NA')

        elif command == 'srog':
            config_log(switch=args.debug, logfile=args.log)
            logging.info(
                "Determine the spacial realtions of genomic (SROG) intervals")
            summary = srog_peak(inbed1=args.bed1,
                                inbed2=args.bed2,
                                outfile=args.output,
                                max_dist=args.max_dist)
            print(summary)

        elif command == 'covary':
            config_log(switch=args.debug, logfile=args.log)
            a_uniq_lst, b_uniq_lst, common_lst = compare_bed(
                args.bed1, args.bed2)
            logging.info("Calculate covariabilities of overlapped regions ...")
            c_corr = bigwig_corr(bed=common_lst,
                                 bw1=args.bw1,
                                 bw2=args.bw2,
                                 outfile=args.output + '_common.tsv',
                                 na_label=args.na_label,
                                 score_type=args.score_type,
                                 exact_scores=args.exact,
                                 keep_NA=args.keepna,
                                 top_x=args.top_X,
                                 min_sig=args.min_signal)
            print(c_corr.T)

            logging.info(
                "Calculate covariabilities of \"%s\" unique regions ..."
                % args.bed1)
            if args.nameA is not None:
                outfile_A = args.output + '_' + args.nameA + '_unique.tsv'
            else:
                outfile_A = args.output + '_bedA_unique.tsv'
            a_corr = bigwig_corr(bed=a_uniq_lst,
                                 bw1=args.bw1,
                                 bw2=args.bw2,
                                 outfile=outfile_A,
                                 na_label=args.na_label,
                                 score_type=args.score_type,
                                 exact_scores=args.exact,
                                 keep_NA=args.keepna,
                                 top_x=args.top_X,
                                 min_sig=args.min_signal)
            print(a_corr.T)

            logging.info(
                "Calculate covariabilities of \"%s\" unique regions ..."
                % args.bed2)
            if args.nameB is not None:
                outfile_B = args.output + '_' + args.nameB + '_unique.tsv'
            else:
                outfile_B = args.output + '_bedB_unique.tsv'
            b_corr = bigwig_corr(bed=b_uniq_lst,
                                 bw1=args.bw1,
                                 bw2=args.bw2,
                                 outfile=outfile_B,
                                 na_label=args.na_label,
                                 score_type=args.score_type,
                                 exact_scores=args.exact,
                                 keep_NA=args.keepna,
                                 top_x=args.top_X,
                                 min_sig=args.min_signal)
            print(b_corr.T)

        elif command == 'cooccur':
            config_log(switch=args.debug, logfile=args.log)
            logging.info(
                "Calculate the co-occurrence of two sets of genomic intervals")
            results = cooccur_peak(inbed1=args.bed1,
                                   inbed2=args.bed2,
                                   name1=args.nameA,
                                   name2=args.nameB,
                                   inbed_bg=args.bed3,
                                   outfile=args.output,
                                   n_cut=args.n_cut,
                                   p_cut=args.p_cut)
            print(results)

        elif command == 'zscore':
            config_log(switch=args.debug, logfile=args.log)
            cal_zscores(args.input, args.output)


if __name__ == '__main__':
    main()
