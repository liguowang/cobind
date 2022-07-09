Stat
============

Description
-------------
Wrapper function. Report basic statistics of genomic intervals, including
 - count
 - total size
 - unique size
 - mean size
 - median size
 - min size
 - max size
 - Standard deviation

and calculate overlapping measurements, including
 - collocation coefficient (C)
 - Jaccard similarity coefficient (J)
 - Sørensen–Dice coefficient (SD)
 - Szymkiewicz–Simpson coefficient (SS)
 - pointwise mutual information (PMI)
 - normalized pointwise mutual information (NPMI)




Usage
-----

:code:`cobind.py stat -h`

::
 
 usage: cobind.py stat [-h] [-b BGSIZE] [-d] input_A.bed input_B.bed
 
 positional arguments:
   input_A.bed           Genomic regions in BED, BED-like or bigBed format. The BED-like
                         format includes: 'bed3', 'bed4', 'bed6', 'bed12', 'bedgraph',
                         'narrowpeak', 'broadpeak', 'gappedpeak'. BED and BED-like format
                         can be plain text, compressed (.gz, .z, .bz, .bz2, .bzip2) or
                         remote (http://, https://, ftp://) files. Do not compress BigBed
                         foramt. BigBed file can also be a remote file.
   input_B.bed           Genomic regions in BED, BED-like or bigBed format. The BED-like
                         format includes: 'bed3', 'bed4', 'bed6', 'bed12', 'bedgraph',
                         'narrowpeak', 'broadpeak', 'gappedpeak'. BED and BED-like format
                         can be plain text, compressed (.gz, .z, .bz, .bz2, .bzip2) or
                         remote (http://, https://, ftp://) files. Do not compress BigBed
                         foramt. BigBed file can also be a remote file.
 
 optional arguments:
   -h, --help            show this help message and exit
   -b BGSIZE, --background BGSIZE
                         The size of the cis-regulatory genomic regions. This is about
                         1.4Gb For the human genome. (default: 1400000000)
   -d, --debug           Print detailed information for debugging.


Example
-------

:code:`cobind.py stat CTCF_ENCFF660GHM.bed RAD21_ENCFF057JFH.bed`

::
  
  2022-07-09 09:44:12 [INFO]  Gathering information for "CTCF_ENCFF660GHM.bed" ...
  2022-07-09 09:44:12 [INFO]  Gathering information for "RAD21_ENCFF057JFH.bed" ...
  A.name                     CTCF_ENCFF660GHM.bed
  A.interval_count                          58684
  A.interval_total_size                  12190325
  A.interval_mean_size                   207.7283
  A.interval_median_size                 240.0000
  A.interval_min_size                          60
  A.interval_max_size                         576
  A.interval_size_SD                      51.5489
  B.name                    RAD21_ENCFF057JFH.bed
  B.interval_count                          33373
  B.interval_total_size                  11381586
  B.interval_mean_size                   341.0417
  B.interval_median_size                 404.0000
  B.interval_min_size                         101
  B.interval_max_size                         553
  B.interval_size_SD                      96.8607
  G.size                          1400000000.0000
  A.size                                 12184840
  Not_A.size                      1387815160.0000
  B.size                                 11130268
  Not_B.size                      1388869732.0000
  A_not_B.size                            7245355
  B_not_A.size                            6190783
  A_and_B.size                            4939485
  A_and_B.exp_size                     96871.8105
  A_or_B.size                            18375623
  Neither_A_nor_B.size            1381624377.0000
  coef.Collocation                         0.4241
  coef.Jaccard                             0.2688
  coef.Dice                                0.4237
  coef.SS                                  0.4438
  A_and_B.PMI                              3.9316
  A_and_B.NPMI                             0.6962
  dtype: object


