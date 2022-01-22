Introduction
============

Genomic features and genome analysis results are generally represented as genomic intervals, for example,
genes, exons, promoters, enhancers, transcription factor binding sites, DNA motifs, CpG islands, 
nucleosomes, heterochromatin, DNA repeats, SNVs, `INDELs <https://en.wikipedia.org/wiki/Indel>`_.

Overlapping two sets of genomic intervals is a routine task in 
genomic analysis. Many tools have been developed to facilitate these operations, such as `Bedtools <https://bedtools.readthedocs.io/en/latest/index.html>`_, `PyRanges <https://github.com/biocore-NTNU/pyranges>`_, and `GenomicRanges <http://www.bioconductor.org/packages/release/bioc/html/GenomicRanges.html>`_.


Overlapping two lists of genomic intervals is NOT like overlapping two lists of gene symbols; the result of overlapping two gene symbols is binary,
while the extent of overlapping between two genomic intervals is a continuous measure. Usually, the *percentage of overlap* is used to quantify the overlapping between two sets of genomic intervals. Although this has been widely used, such an oversimplified measurement has several serious drawbacks. 

- Problem 1: Subjective. We have to *arbitrarily* define a threshold in order to call if two genomic intervals are overlapped or not. One nucleotide overlap? Ten nucleotide overlap? or 10% overlap?

- Problem 2: Non-symmetric. Assuming we identified 1000 TFBSs for protein A and 10,000 TFBSs for protein B, and we found 500 sites were overlapped with at least one nucleotide. We found 50% (500/1000) of A's binding sites were overlapped with B, while only 5% (500/10000) of B's binding sites were overlapped with A. And the situation becomes even worse when the sizes of two overlapped genomic intervals are significantly different. 

- Problem 3: Incomparable. Suppose we did another analysis and found only 20% of A's binding sites were overlapped with **protein C**. Can we reach the conclusion 
  that A's binding sites overlap better with B than C? The answer is "no" because B and C might have different numbers of TFBS. The answer is still "no" even
  B and C have the same number of TFBS, because the length of TFBS might be significantly different between B and C. 

These drawbacks make the direct comparison (or meta-analysis) of genomic overlapping analysis results problematic. `Cobind <https://cobind.readthedocs.io/en/latest/>`_ is designed to address these problems and provides objective, rigorous quantification of a genomic region overlapping.

