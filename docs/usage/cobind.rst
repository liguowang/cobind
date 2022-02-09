Overview
=========

Subcommands description
-----------------------
`cobind <https://cobind.readthedocs.io/en/latest/index.html>`_ is a python package designed to quantify the "overlapping" of genomic intervals.

.. list-table:: **subcommands provided by cobind**
   :widths: 15,50
   :header-rows: 1

   * - *Subcommand*
     - Description
   * - `overlap <https://cobind.readthedocs.io/en/latest/usage/overlap.html>`_
     - Calculate the `overlap coefficient (O) <https://cobind.readthedocs.io/en/latest/definition.html#overlap-coefficient-o>`_.
   * - `jaccard <https://cobind.readthedocs.io/en/latest/usage/jaccard.html>`_
     - Calculate the `Jaccard similarity coefficient (J) <https://cobind.readthedocs.io/en/latest/definition.html#jaccard-coefficient-j>`_.
   * - `dice <https://cobind.readthedocs.io/en/latest/usage/SD.html>`_
     - Calculate the `Sørensen–Dice coefficient (SD) <https://cobind.readthedocs.io/en/latest/definition.html#sorensendice-coefficient-sd>`_.
   * - `simpson <https://cobind.readthedocs.io/en/latest/usage/SS.html>`_
     - Calculate the `Szymkiewicz–Simpson coefficient (SS) <https://cobind.readthedocs.io/en/latest/definition.html#szymkiewiczsimpson-coefficient-ss>`_.
   * - `pmi <https://cobind.readthedocs.io/en/latest/usage/PMI.html>`_
     - Calculate the `pointwise mutual information (PMI) <https://cobind.readthedocs.io/en/latest/definition.html#pointwise-mutual-information-pmi>`_.
   * - `npmi <https://cobind.readthedocs.io/en/latest/usage/NPMI.html>`_
     - Calculate the `normalized pointwise mutual information (NPMI) <https://cobind.readthedocs.io/en/latest/definition.html#normalized-pointwise-mutual-information-npmi>`_.
   * - `cooccur <https://cobind.readthedocs.io/en/latest/usage/cooccur.html>`_
     - Evaluate if two sets of genomic regions are significantly overlapped.
   * - `covary <https://cobind.readthedocs.io/en/latest/usage/covary.html>`_
     - Calculate the covariance of binding intensities between two sets of genomic intervals.
   * - `srog <https://cobind.readthedocs.io/en/latest/usage/SROG.html>`_
     - Report the code of `Spatial Relation Of Genomic (SROG) <https://cobind.readthedocs.io/en/latest/definition.html#spacial-relations-of-genomic-regions-srog>`_ regions.
   * - `stat <https://cobind.readthedocs.io/en/latest/usage/stat.html>`_
     - Wrapper function. Calculate *O*, *J*, *SD*, *SS*, *PMI*, and *NPMI*.



Usage
-----

Print out all the avaiable commands

:code:`cobind.py -h`

::
 
 usage: cobind.py [-h] [-v] {overlap,jaccard,dice,simpson,pmi,npmi,cooccur,covary,srog,stat} ...
 
 **cobind: collocation analyses of genomic intervals**
 
 positional arguments:
   {overlap,jaccard,dice,simpson,pmi,npmi,cooccur,covary,srog,stat}
                         Sub-command description:
     overlap             Calculate the overlapping coefficient (O) between two sets of genomic
                         regions. O = |A and B| / (|A|*|B|)**0.5
     jaccard             Calculate the Jaccard similarity coefficient (J) between two sets of genomic
                         regions. J = |A and B| / |A or B|
     dice                Calculate the Sørensen–Dice coefficient (SD) between two sets of genomic
                         regions. SD = 2*|A and B| / (|A| + |B|)
     simpson             Calculate the Szymkiewicz–Simpson coefficient (SS) between two sets of
                         genomic regions. SS = |A and B| / min(|A|, |B|)
     pmi                 Calculate the pointwise mutual information (PMI) between two sets of genomic
                         regions. PMI = log(p(|A and B|)) - log(p(|A|)) - log(p(|B|))
     npmi                Calculate the normalized pointwise mutual information (NPMI) between two sets
                         of genomic regions. NPMI = log(p(|A|)*p(|B|)) / log(p(|A and B|)) - 1
     cooccur             Evaluate if two sets of genomic regions are significantly overlapped in given
                         background regions.
     covary              Calculate the covariance (Pearson, Spearman and Kendall coefficients) of
                         binding intensities between two sets of genomic regions.
     srog                Report the code of Spatial Relation Of Genomic (SROG) regions. SROG codes
                         include 'disjoint','touch','equal','overlap','contain','within'.
     stat                Wrapper function. Report basic statistics of genomic regions, and calculate
                         overlapping measurements, including "O", "J", "SD", "SS", "PMI", "NPMI",
                         without bootstrap resampling or generating peakwise measurements.
 
 optional arguments:
   -h, --help            show this help message and exit
   -v, --version         show program's version number and exit

