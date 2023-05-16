Definitions
============


Symbols definitions
-------------------
Say, we have two sets of genomic intervals **A** and **B**, and the genomic background is **G**.
In Figure 1 below, both A and B contain only 1 genomic region for the purpose of clarity, but the definitions
are still applicable if **A** and **B** contain many intervals. 

Symbols are defined as:

.. image:: _static/set_symbols.jpg
  :width: 600
  :alt: Figure_1


\|A\|
  The `cardinality <https://en.wikipedia.org/wiki/Cardinality>`_ of **A** (i.e., all the **non-redundant** bases covered by **A**). For example, if A contains two genomic intervals: "chr1 0 10", "chr1 5 15", then \|A\| = 15. 
\|B\|
  The `cardinality <https://en.wikipedia.org/wiki/Cardinality>`_ of **B** (i.e., all the **non-redundant** bases covered by **B**).
\|G\|
  The genomic background. Depending on the context, this can be *the whole genome*, *all the cis-regulatory elements*, *all the promoters*, *all the TF binding sites* in the genome, etc. **A** and **B** must be the subsets of **G**. 
\|A ∪ B\|
  Union of A and B (i.e., bases covered by A or B).
\|A ∩ B\|
  Intersection of A and B (i.e., bases covered by A and B simultaneously). This is commonly used to measure the *collocation* of A and B.
\|A − B\|
  Difference (A not B) (i.e., bases covered by only A but not B).
\|B − A\|
  Difference (B not A) (i.e., bases covered by only B but not A).
\|A ∪ B\|^𝐶
  Complement of \|A ∪ B\| (i.e., bases NOT covered by A or B).


Spacial Relations of Genomic regions (SROG)
-------------------------------------------

There are six different spacial relations between two genomic regions (A and B). As illustrated below:

.. image:: _static/srog.jpg
  :width: 700
  :alt: Alternative text


Collocation coefficient (C)
---------------------------
The collocation coefficient between A and B is calculated as the ratio between \|A ∩ B\| and the *geometric mean of \|A\| and \|B\|*.
C(A,B) is a value between [0, 1], with 0 indicating 'no overlap', and 1 indicating '100% overlap' (i.e., A and B are identical). C(A, B) is defined as 0 when \|A\| = 0 or \|B\| = 0, or  \|A\| = \|B\| = 0.

.. image:: _static/ov_coef_1.jpg
  :width: 250
  :alt: Alternative text

.. image:: _static/ov_coef_3.jpg
  :width: 200
  :alt: Alternative text

Overall collocation coefficient
  The collocation coefficient between two **sets** of genomic regions. For example, you can use the *overall collocation coefficient* to measure the cobindability of two transcription factors. 

peakwise collocation coefficient
  The collocation coefficient between **two** genomic intervals (A protein-bound genomic region is called "peak" in `ChIP-seq <https://en.wikipedia.org/wiki/ChIP_sequencing>`_ experiment). 



Jaccard coefficient (J)
-------------------------
The `Jaccard similarity coefficient <https://en.wikipedia.org/wiki/Jaccard_index>`_, also known as the Jaccard index. It is the ratio between **intersection** and **union**. J(A, B) is defined as 0 when \|A\| = 0 or \|B\| = 0, or  \|A\| = \|B\| = 0.


.. image:: _static/jaccard_1.jpg
  :width: 400
  :alt: Alternative text

.. image:: _static/jaccard_2.jpg
  :width: 180
  :alt: Alternative text


The Jaccard distance *Dj* is calculated as:

.. image:: _static/jaccard_3.jpg
  :width: 450
  :alt: Alternative text


Similar to O(A,B), we have an **overall Jaccard coefficient** and **peakwise Jaccard coefficient**.

.. note::
   The Jaccard coefficient implemented here is slightly different from `BEDTools <https://bedtools.readthedocs.io/en/latest/content/tools/jaccard.html>`_ :code:`jaccard` function.
   When calculating the union, BEDTools only use the intervals that are overlapped with each other, while we use all the intervals.

overall Jaccard coefficient
  The Jaccard coefficient between two **sets** of genomic regions. 
peakwise Jaccard coefficient
  The Jaccard coefficient between **two** genomic intervals.


Sørensen–Dice coefficient (SD)
------------------------------
`Sørensen–Dice coefficient <https://en.wikipedia.org/wiki/S%C3%B8rensen%E2%80%93Dice_coefficient>`_,  also called *Sørensen–Dice index*, *Sørensen index* or *Dice's coefficient*. SD(A, B) is defined as 0 when \|A\| = 0 or \|B\| = 0, or  \|A\| = \|B\| = 0.

.. image:: _static/SD_1.jpg
  :width: 200
  :alt: Alternative text

.. image:: _static/SD_2.jpg
  :width: 180
  :alt: Alternative text

Jaccard coefficient (J) can be converted into Sørensen–Dice coefficient (SD) and vice versa:

*J = SD/(2-SD)* and *SD = 2J/(1+J)*


Szymkiewicz–Simpson coefficient (SS)
-------------------------------------
`Szymkiewicz–Simpson coefficient <https://en.wikipedia.org/wiki/Overlap_coefficient>`_ is defined as the size of the intersection divided by the smaller of the size of the two sets.

.. image:: _static/SS.jpg
  :width: 250
  :alt: Alternative text

.. image:: _static/SS_bound.jpg
  :width: 180
  :alt: Alternative text


Pointwise mutual information (PMI)
----------------------------------
`Pointwise mutual information (PMI) <https://en.wikipedia.org/wiki/Pointwise_mutual_information>`_ is one of the standard association measures in collocation analysis. 
It measures how much the observed overlaps differ from what we would expect them to be. Assume A and B represent two sets of genomic regions bound by `transcription factors <https://en.wikipedia.org/wiki/Transcription_factor>`_ A and B; respectively, PMI measures if A and B bind together or separately.


PMI is calculated as:

.. image:: _static/pmi.jpg
  :width: 300
  :alt: Alternative text

where 

.. image:: _static/p.jpg
  :width: 300
  :alt: Alternative text


PMI = 0
  Indicates that A and B are independent.
PMI > 0
  Indicates that the overlapping between A and B is in a frequency *higher* than what we would expect if A and B are independent (i.e, A and B tend to bind together). 
PMI < 0
  Indicates that the overlapping between A and B is in frequency *lower* than what we would expect if A and B are independent. (i.e., A and B tend to bind separately). 

Note, PMI has no boundaries:

.. image:: _static/pmi_bound.jpg
  :width: 500
  :alt: Alternative text


Normalized pointwise mutual information (NPMI)
----------------------------------------------
Normalized pointwise mutual information (NPMI) is calculated as:

.. image:: _static/npmi.jpg
  :width: 650
  :alt: Alternative text

Note, after normalization, NPMI is confined to [-1, 1]:

.. image:: _static/npmi_bound.jpg
  :width: 250
  :alt: Alternative text


Which metric to use? 
---------------------

Based on our evaluation, the **Collocation coefficient (C)** and **NPMI** are the best two metrics one can use to quantify the overlap (collocation)
between two sets of genomic intervals.

`Metric evaluation <https://cobind.readthedocs.io/en/latest/comparison.html>`_



