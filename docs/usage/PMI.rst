Pointwise mutual information (PMI)
====================================

Description
-------------

Calculate the `Pointwise mutual information (PMI) <https://en.wikipedia.org/wiki/Pointwise_mutual_information>`_ [#f1]_ between two sets of genomic regions. 

.. image:: ../_static/pmi.jpg
  :width: 300
  :alt: Alternative text

.. image:: ../_static/pmi_bound.jpg
  :width: 500
  :alt: Alternative text

where 

.. image:: ../_static/p.jpg
  :width: 300
  :alt: Alternative text



Command (getting help)
----------------------

:code:`cobind.py pmi -h`

::

 usage: cobind.py pmi [-h] [-n ITER] [-f SUBSAMPLE] [-b BGSIZE] [-o] [-d] input_A.bed input_B.bed
 
  positional arguments:
    input_A.bed           Genomic regions in BED, BED-like or bigBed format. The BED-like format
                          includes: 'bed3', 'bed4', 'bed6', 'bed12', 'bedgraph', 'narrowpeak',
                          'broadpeak', 'gappedpeak'. BED and BED-like format can be plain text,
                          compressed (.gz, .z, .bz, .bz2, .bzip2) or remote (http://, https://,
                          ftp://) files. Do not compress BigBed foramt. BigBed file can also be a
                          remote file.
    input_B.bed           Genomic regions in BED, BED-like or bigBed format. The BED-like format
                          includes: 'bed3', 'bed4', 'bed6', 'bed12', 'bedgraph', 'narrowpeak',
                          'broadpeak', 'gappedpeak'. BED and BED-like format can be plain text,
                          compressed (.gz, .z, .bz, .bz2, .bzip2) or remote (http://, https://,
                          ftp://) files. Do not compress BigBed foramt. BigBed file can also be a
                          remote file.
  
  optional arguments:
    -h, --help            show this help message and exit
    -n ITER, --ndraws ITER
                          Times of resampling to estimate confidence intervals. Set to '0' to turn
                          off resampling.(default: 20)
    -f SUBSAMPLE, --fraction SUBSAMPLE
                          Resampling fraction. (default: 0.75)
    -b BGSIZE, --background BGSIZE
                          The size of the cis-regulatory genomic regions. This is about 1.4Gb For the
                          human genome. (default: 1400000000)
    -o, --save            If set, will save peak-wise coefficients to files
                          ("input_A_peakwise_scores.tsv" and "input_B_peakwise_scores.tsv").
    -d, --debug           Print detailed information for debugging.


Command (example)
-----------------

Calculate the **overall** `PMI <https://en.wikipedia.org/wiki/Pointwise_mutual_information>`_ and **peak-wise** `PMI <https://en.wikipedia.org/wiki/Pointwise_mutual_information>`_ between CTCF binding sites and RAD21 binding sites.

:code:`python3 ../bin/cobind.py pmi CTCF_ENCFF660GHM.bed RAD21_ENCFF057JFH.bed --save`

The overall `PMI <https://en.wikipedia.org/wiki/Pointwise_mutual_information>`_ between *CTCF_ENCFF660GHM.bed* and *RAD21_ENCFF057JFH.bed* was printed to screen

::

 2022-01-16 09:01:34 [INFO]  Calculate the pointwise mutual information (PMI) ...
 A.name               CTCF_ENCFF660GHM.bed
 B.name              RAD21_ENCFF057JFH.bed
 A.interval_count                    58684
 B.interval_count                    33373
 A.size                           12184840
 B.size                           11130268
 A_or_B.size                      18375623
 A_and_B.size                      4939485
 Coef                               3.9316
 Coef(expected)                     0.0000
 Coef(95% CI)              [3.9230,3.9343]
 dtype: object
 2022-01-16 09:02:02 [INFO]  Read and union BED file: "CTCF_ENCFF660GHM.bed"
 2022-01-16 09:02:03 [INFO]  Unioned regions of "CTCF_ENCFF660GHM.bed" : 58584
 2022-01-16 09:02:03 [INFO]  Read and union BED file: "RAD21_ENCFF057JFH.bed"
 2022-01-16 09:02:03 [INFO]  Unioned regions of "RAD21_ENCFF057JFH.bed" : 31955
 2022-01-16 09:02:03 [INFO]  Build interval tree for unioned BED file: "CTCF_ENCFF660GHM.bed"
 2022-01-16 09:02:03 [INFO]  Build interval tree for unioned BED file: "RAD21_ENCFF057JFH.bed"
 2022-01-16 09:02:03 [INFO]  Calculate the overlap coefficient of each genomic region in CTCF_ENCFF660GHM.bed ...
 2022-01-16 09:02:06 [INFO]  Save peakwise scores to CTCF_ENCFF660GHM.bed_peakwise_scores.tsv ...
 2022-01-16 09:02:06 [INFO]  Calculate the overlap coefficient of each genomic region in RAD21_ENCFF057JFH.bed ...
 2022-01-16 09:02:07 [INFO]  Save peakwise scores to RAD21_ENCFF057JFH.bed_peakwise_scores.tsv ...

If :code:`--save` was specified, the peakwise `PMI <https://en.wikipedia.org/wiki/Pointwise_mutual_information>`_ were saved to *CTCF_ENCFF660GHM.bed_peakwise_scores.tsv* and *RAD21_ENCFF057JFH.bed_peakwise_scores.tsv*, respectively.
::

 $ head -5 CTCF_ENCFF660GHM.bed_peakwise_scores.tsv
  
 chrom start end A.size  B.size  A∩B A∪B B.list  Score
 chr12 108043  108283  240 404 240 404 chr12:107919-108323 15.058323195606475
 chr12 153232  153470  238 222 222 238 chr12:153236-153458 15.58746739989615
 chr12 177749  177989  240 NA  NA  NA  NA  NA
 chr12 189165  189405  240 404 240 404 chr12:189072-189476 15.058323195606475

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
  The peakwise `PMI <https://en.wikipedia.org/wiki/Pointwise_mutual_information>`_.


.. [#f1] The natural log was used when calculating PMI.
