Introduction
============

Most genomic features and genome analysis results are represented as genomic intervals (e.g., 
genes, exons, promoters, enhancers, transcription factor binding sites, DNA motifs, CpG islands, 
nucleosomes, heterochromatin, Alu repeats). Single nucleotide variants (SNVs), and `INDELs <https://en.wikipedia.org/wiki/Indel>`_ can also represent as
genomic intervals.

Overlapping two sets of genomic features (represented as genomic intervals) is one of the most common tasks in 
genomic analysis. For example, you may want to overlap two sets of transcription factor binding sites (TFBSs)
identified from `ChIP-seq <https://en.wikipedia.org/wiki/ChIP_sequencing>`_ experiments or overlap TFBSs with
known genomic features (genes, exons, promoters, enhancers, CpG islands) to annotate them. 
Many tools have been developed to facilitate these operations, such as `Bedtools <https://bedtools.readthedocs.io/en/latest/index.html>`_, `PyRanges <https://github.com/biocore-NTNU/pyranges>`_, and `GenomicRanges <http://www.bioconductor.org/packages/release/bioc/html/GenomicRanges.html>`_.



Overlap two lists of genomic intervals is NOT like overlap two lists of gene symbols; the result of overlapping two gene symbols is binary,
while the extent of overlapping between two genomic intervals is continuous. Typically, the *percentage of overlap* is used to quantify the overlapping between two sets of genomic intervals. Although this has been widely used, such a oversimplified approach has several serious drawbacks. Without loss of generality, let us assume we identified 1000 TFBSs for protein A and 10,000 TFBSs for protein B, and we found 500 sites were overlapped with at least one nucleotide.

- Problem 1: Arbitrary. We have to *arbitrarily* define a threshold before doing the overlapping. One nucleotide overlap? Ten nucleotide overlap?
  or 10% overlap? Is this criterion reciprocal? A reciprocal cutoff is also problematic when one genomic interval is significantly larger than another. 
- Problem 2: Non-symmetric.  50% (500/1000) of A's binding sites were overlapped with B, while only 5% (500/10000) of B's binding sites were overlapped with A. 
  Which *percentage of overlap* to report?
- Problem 3: Incomparable. Suppose we did another analysis and found only 20% of A's binding sites were overlapped with **protein C**. Can we reach the conclusion 
  that A's binding sites overlap better with B than C? The answer is "no" because B and C might have different number of TFBS. The answer is still "no" even
  B and C have the same number of TFBS, because the length of TFBS might be significantly different between B and C. 

These problems make the direct comparison of genomic overlapping analysis results difficult even impossible. Cobind is designed to address these problems and provides rigorous quantification of a genomic region overlapping.

