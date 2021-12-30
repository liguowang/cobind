overlap
========

Description
-------------

Calculate the overlap coefficient between two sets of genomic regions. 

Command (getting help)
----------------------

:code:`cobind.py overlap -h`

::

 usage: cobind.py overlap [-h] [-n ITER] [-s SUBSIZE] [-b BGSIZE] [-o] [-d] input_A.bed input_B.bed
 
 positional arguments:
   input_A.bed           Genomic regions in BED, BED-like or bigBed format. The BED-like format includes: 'bed3', 'bed4', 'bed6', 'bed12', 'bedgraph', 'narrowpeak', 'broadpeak', 'gappedpeak'. BED and BED-like format can be plain text, compressed (.gz, .z,
                        .bz, .bz2, .bzip2) or remote (http://, https://, ftp://) files. Do not compress BigBed foramt. BigBed file can also be a remote file.
   input_B.bed           Genomic regions in BED, BED-like or bigBed format. The BED-like format includes: 'bed3', 'bed4', 'bed6', 'bed12', 'bedgraph', 'narrowpeak', 'broadpeak', 'gappedpeak'. BED and BED-like format can be plain text, compressed (.gz, .z,
                         .bz, .bz2, .bzip2) or remote (http://, https://, ftp://) files. Do not compress BigBed foramt. BigBed file can also be a remote file.
 
 optional arguments:
  -h, --help            show this help message and exit
  -n ITER, --ndraws ITER
                        Times of resampling to estimate confidence intervals. Set to '0' to turn off resampling.(default: 50)
  -s SUBSIZE, --size SUBSIZE
                        Size of the subset during resampling. If the original BED file contains 10000 regions, '--size = 0.85' means 8500 regions will be resampled. (default: 0.85)
  -b BGSIZE, --background BGSIZE
                        The size of the cis-regulatory genomic regions. This is about 1.4Gb For the human genome. (default: 1400000000)
  -o, --save            If set, peak-wise overlapping coefficients will save to files ("input_A_ovcoef.tsv" and "input_B_ovcoef.tsv").
  -d, --debug           Print detailed information for debugging.


Command (example)
-----------------

Calculate the **overall** overlap coefficient and **peak-wise** overlap coefficients between CTCF binding sits and RAD21 binding sites.

:code:`python3 ../bin/cobind.py overlap CTCF_ENCFF660GHM.bed RAD21_ENCFF057JFH.bed --save`


The overall overlapping coefficient between *CTCF_ENCFF660GHM.bed* and *RAD21_ENCFF057JFH.bed* was 
printed to screen:

- coef_obs : observed overlapc oefficient. Please note, *0 <= coef_obs <= 1* with "0" indicating no overlap and "1" indicating complete overlap (i.e., two BED files are identical!)
- coef_exp : expected overlap coefficient.
- coef_ratio : ratio between *coef_obs* and *coef_exp*.
- coef_ratio_low : lower bound of 95% confidence interval of *coef_ratio*.
- coef_ratio_high : upper bound of 95% confidence interval of *coef_ratio*.


::

 2021-12-29 08:09:03 [INFO]  Calculate overall overlapping coefficient ...
 coef_obs            0.424149
 coef_exp            0.008318
 coef_ratio         50.989911
 coef_ratio_low     50.691262
 coef_ratio_high    51.222445
 dtype: float64
 2021-12-29 08:10:14 [INFO]  Calculate peak-wise overlapping coefficient ...
 2021-12-29 08:10:14 [INFO]  Read and union BED file: "CTCF_ENCFF660GHM.bed"
 2021-12-29 08:10:15 [INFO]  Unioned regions of "CTCF_ENCFF660GHM.bed" : 58584
 2021-12-29 08:10:15 [INFO]  Read and union BED file: "RAD21_ENCFF057JFH.bed"
 2021-12-29 08:10:15 [INFO]  Unioned regions of "RAD21_ENCFF057JFH.bed" : 31955
 2021-12-29 08:10:15 [INFO]  Build interval tree for unioned BED file: "CTCF_ENCFF660GHM.bed"
 2021-12-29 08:10:15 [INFO]  Build interval tree for unioned BED file: "RAD21_ENCFF057JFH.bed"
 2021-12-29 08:10:15 [INFO]  Calculate the overlap coefficient of each genomic region in CTCF_ENCFF660GHM.bed ...
 2021-12-29 08:10:17 [INFO]  Calculate the overlap coefficient of each genomic region in RAD21_ENCFF057JFH.bed ...
 

Check peak-wise overlap coefficients. Using the *CTCF_ENCFF660GHM_ovcoef.tsv* file as an example,
there are a total of 8 columns in this file.

- column 1 to 3: genomic coordinate of CTCF peak.
- column 4: The number of RAD21 peaks overlapped with this peak.
- column 5: The number of bases overlapped with this peak.
- column 6: The ratio between "number of bases overlapped" and "total bases of the peak". 
- column 7: overlap coefficient.
- column 8: RAD21 peak(s) that are overlapped with this peak. Multiple peaks will be separated by ",".

::

 $head -5 CTCF_ENCFF660GHM_ovcoef.tsv
 
 #chrom  start end ov_peaks_n  ov_bases_n  ov_bases_frac ov_coef ov_peaks_list
 chr12 108043  108283  1 240 1.0 0.770752493308062 chr12:107919-108323
 chr12 153232  153470  1 222 0.9327731092436975  0.965801796044974 chr12:153236-153458
 chr12 177749  177989  0 NA  NA  NA  NA
 chr12 189165  189405  1 240 1.0 0.770752493308062 chr12:189072-189476
  
 
 $head -5 RAD21_ENCFF057JFH_ovcoef.tsv
 
 #chrom  start end ov_peaks_n  ov_bases_n  ov_bases_frac ov_coef ov_peaks_list
 chr22 16146910  16147314  1 240 0.594059405940594 0.770752493308062 chr22:16147004-16147244
 chr22 16391500  16391904  1 90  0.22277227722772278 0.47198758164566446 chr22:16391656-16391746
 chr22 16441033  16441408  1 240 0.64  0.8 chr22:16441063-16441303
 chr22 16599668  16600072  0 NA  NA  NA  NA

