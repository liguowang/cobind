Introduction
============

Introductionn
-------------

Collocated genomic intervals indicate biological association. Therefore, overlapping analysis of genomic intervals has been widely used to QC, integrate, and impute the function of genomic intervals. 

The conventional approach to measure the "overlap between genomic intervals" involves arbitrary thresholds to decide the total number of overlapped genomic regions, which leads to biased, non-reproducible, and incomparable results. Specifically,

 * The result derived from this *threshold-and-count* approach is non-reproducible and incomparable, as different thresholds produce different results.
 * The overlapping between two intervals is a continuous variable, whereas the thresholded approach reduces it into a binary variable. Casting the one-dimensional intervals as zero-dimensional points loses the information and sensitivity needed to accurately evaluate the collocation strength.
 * The absolute or relative counts is biased by the size and the total number of intervals.

Cobind provides six threshold-free metrics to quantify the *strength of genomic overlapping* rigorously. 

 * `the Collocation coefficient (C) <https://cobind.readthedocs.io/en/latest/definition.html#collocation-coefficient-c>`_
 * `the Jaccard coefficient (J) <https://cobind.readthedocs.io/en/latest/definition.html#jaccard-coefficient-j>`_
 * `the Sørensen–Dice coefficient (SD) <https://cobind.readthedocs.io/en/latest/definition.html#sorensendice-coefficient-sd>`_
 * `the Szymkiewicz–Simpson coefficient (SS) <https://cobind.readthedocs.io/en/latest/definition.html#szymkiewiczsimpson-coefficient-ss>`_
 * `the Pointwise Mutual Information (PMI) <https://cobind.readthedocs.io/en/latest/definition.html#pointwise-mutual-information-pmi>`_
 * `the Normalized Pointwise Mutual Information (NPMI) <https://cobind.readthedocs.io/en/latest/definition.html#normalized-pointwise-mutual-information-npmi>`_
