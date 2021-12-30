.. role:: raw-math(raw)
    :format: latex html

Input file and data format
===========================

BED format
----------
`BED (Browser Extensible Data) <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>`_ 
format is commonly used to describe genomic intervals. BED files used by **cobind** must 
have at least the first 3 columns::

 # BED3 format (chrom, start, end)
 chr1    629149    629391
 chr1    629720    630165
 chr1    631404    631758
 ...
 
 # BED4 format (chrom, start, end, name)
 chr1    629149    629391   region_1
 chr1    629720    630165   region_2
 chr1    631404    631758   region_3
 ...
 
 # BED6 format (chrom, start, end, name, score, strand)
 chr1    629149  629391 region_1    0    +
 chr1    629720  630165 region_2    0    +
 chr1    631404  631758 region_3    0    -
 ...

BED-like format
---------------

- `bedgraph <https://genome.ucsc.edu/goldenPath/help/bedgraph.html>`_
- ENCODE `narrowpeak <https://genome.ucsc.edu/FAQ/FAQformat.html#format12>`_
- ENCODE `broadpeak <https://genome.ucsc.edu/FAQ/FAQformat.html#format13>`_
- ENCODE `gappedpeak <https://genome.ucsc.edu/FAQ/FAQformat.html#format14>`_


bigBed
------
`bigBed <https://genome.ucsc.edu/goldenPath/help/bigBed.html>`_ is an an indexed binary format of BED file. `UCSC's <http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/>`_  :code:`bedToBigBed` and :code:`bigBedToBed` commands can be used to convert BED file into bigBed file or *vice versa*.


bigWig
------
The `bigWig <https://genome.ucsc.edu/goldenpath/help/bigWig.html>`_ format is an indexed binary format of `wiggle <https://genome.ucsc.edu/goldenpath/help/wiggle.html>`_ file, which is widely used to represent genomic signals. `UCSC's <http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/>`_  :code:`wigToBigWig` and :code:`bigWigToWig` commands can be used to convert wiggle file into bigWig file or *vice versa*.

