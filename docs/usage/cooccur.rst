Co-occurrence
==============

Description
-------------

Evaluate if two sets of genomic regions are significantly overlapped in given background regions.

.. image:: ../_static/ov_coef_1.jpg
  :width: 250
  :alt: Alternative text

.. image:: ../_static/ov_coef_3.jpg
  :width: 200
  :alt: Alternative text

Command (getting help)
----------------------

:code:`cobind.py cooccur -h`

::

 usage: cobind.py cooccur [-h] [--ncut N_CUT] [--pcut P_CUT] [-d] [--dist MAX_DIST] input_A.bed input_B.bed background.bed output.tsv
 
 positional arguments:
   input_A.bed      Genomic regions in BED, BED-like or bigBed format. The BED-like format includes: 'bed3', 'bed4', 'bed6', 'bed12',
                    'bedgraph', 'narrowpeak', 'broadpeak', 'gappedpeak'. BED and BED-like format can be plain text, compressed (.gz, .z,
                    .bz, .bz2, .bzip2) or remote (http://, https://, ftp://) files. Do not compress BigBed foramt. BigBed file can also
                    be a remote file.
   input_B.bed      Genomic regions in BED, BED-like or bigBed format. The BED-like format includes: 'bed3', 'bed4', 'bed6', 'bed12',
                    'bedgraph', 'narrowpeak', 'broadpeak', 'gappedpeak'. BED and BED-like format can be plain text, compressed (.gz, .z,
                    .bz, .bz2, .bzip2) or remote (http://, https://, ftp://) files. Do not compress BigBed foramt. BigBed file can also
                    be a remote file.
   background.bed   Genomic regions as the background (e.g., all promoters, all enhancers).
   output.tsv       For each genomic region in the background BED file, add another column indicating if this region is "input_A
                    specific", "input_B specific", "co-occur" or "neither".
 
 optional arguments:
   -h, --help       show this help message and exit
   --ncut N_CUT     The minimum overlap size. default: 1)
   --pcut P_CUT     The minimum overlap percentage. default: 0.000000)
   -d, --debug      Print detailed information for debugging.
   --dist MAX_DIST  When intervals are disjoint, find the closest up- and down-stream intervals that are no further than `max_dist` away.
                    default: 250000000)



Command (example)
-----------------

Calculate the **overall** overlap coefficient and **peak-wise** overlap coefficients between CTCF binding sites and RAD21 binding sites.

:code:`python3 ../bin/cobind.py overlap CTCF_ENCFF660GHM.bed RAD21_ENCFF057JFH.bed --save`

The overall overlapping coefficient between *CTCF_ENCFF660GHM.bed* and *RAD21_ENCFF057JFH.bed* was printed to screen

::

 2022-01-16 07:47:15 [INFO]  Calculate overlapping coefficient (overall) ...
 A.name               CTCF_ENCFF660GHM.bed
 B.name              RAD21_ENCFF057JFH.bed
 A.interval_count                    58684
 B.interval_count                    33373
 A.size                           12184840
 B.size                           11130268
 A_or_B.size                      18375623
 A_and_B.size                      4939485
 Coef                               0.4241
 Coef(expected)                     0.0083
 Coef(95% CI)              [0.4223,0.4269]
 dtype: object
 2022-01-16 07:47:42 [INFO]  Calculate overlapping coefficient (peak-wise) ...
 2022-01-16 07:47:42 [INFO]  Read and union BED file: "CTCF_ENCFF660GHM.bed"
 2022-01-16 07:47:42 [INFO]  Unioned regions of "CTCF_ENCFF660GHM.bed" : 58584
 2022-01-16 07:47:42 [INFO]  Read and union BED file: "RAD21_ENCFF057JFH.bed"
 2022-01-16 07:47:43 [INFO]  Unioned regions of "RAD21_ENCFF057JFH.bed" : 31955
 2022-01-16 07:47:43 [INFO]  Build interval tree for unioned BED file: "CTCF_ENCFF660GHM.bed"
 2022-01-16 07:47:43 [INFO]  Build interval tree for unioned BED file: "RAD21_ENCFF057JFH.bed"
 2022-01-16 07:47:43 [INFO]  Calculate the overlap coefficient of each genomic region in CTCF_ENCFF660GHM.bed ...
 2022-01-16 07:47:45 [INFO]  Save peakwise scores to CTCF_ENCFF660GHM.bed_peakwise_scores.tsv ...
 2022-01-16 07:47:45 [INFO]  Calculate the overlap coefficient of each genomic region in RAD21_ENCFF057JFH.bed ...
 2022-01-16 07:47:46 [INFO]  Save peakwise scores to RAD21_ENCFF057JFH.bed_peakwise_scores.tsv ...
 


If :code:`--save` was specified, the peakwise overlap coefficients were saved to *CTCF_ENCFF660GHM.bed_peakwise_scores.tsv* and *RAD21_ENCFF057JFH.bed_peakwise_scores.tsv*, respectively.
::

 $ head -5 CTCF_ENCFF660GHM.bed_peakwise_scores.tsv
  
 chrom start end A.size  B.size  A∩B A∪B B.list  Score
 chr12 108043  108283  240 404 240 404 chr12:107919-108323 0.770752493308062
 chr12 153232  153470  238 222 222 238 chr12:153236-153458 0.965801796044974
 chr12 177749  177989  240 NA  NA  NA  NA  NA
 chr12 189165  189405  240 404 240 404 chr12:189072-189476 0.770752493308062

column 1 to 3
  The genomic coordinate of CTCF peak.
column 4 (A.size)
  The size of CTCF peak.
column 5 (B.size)
  The size (cardinality) of RAD21 peak(s) that were overlapped with this CTCF peak.
column 6 (A∩B)
  The size (cardinality) of intersection.
column 7 (A∪B)
  The size (cardinality) of union.
column 8 (B.list)
  List of RAD21 peak(s) that are overlapped with this peak. Multiple peaks will be separated by ",".
column 9 (Score)
  The peakwise overlap coefficient.


.. [#f1] Do not confuse with `Szymkiewicz–Simpson coefficient <https://en.wikipedia.org/wiki/Overlap_coefficient>`_, which is called "overlap coefficent" in Wikipedia, but was named as the "SS coefficient" in our cobind package.