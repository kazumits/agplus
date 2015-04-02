agplus: a rapid and flexible tool for aggregation plots
===========================================================

What's this for?
----------------

**agplus**, is a simple command-line tool, enables rapid and flexible production of text tables tailored for aggregation plots from which users can easily design multiple groups based on user-definitions such as regulatory regions or transcription initiation sites.


Installation
------------

Ruby and BEDTools are required to run **agplus**.

* [Ruby](https://www.ruby-lang.org) (version >2.0)
* [BEDTools](http://bedtools.readthedocs.org)

Place the **agplus** file and the supporting shell scripts anywhere you like. You will need the below software to use our supporting shell scripts (**bam2bwshifted**, **agpdraw-line** and **assignExprGroupPer10**) and to follow our Tutorial at the last section in this README.

* wigToBigWig: required by **bam2bwshifted**
* bigWigToWig: required in the Tutorial
* [SAMtools](http://samtools.sourceforge.net): required by **bam2bwshifted**
* [R](http://www.r-project.org): required by **agpdraw-line** and **assignExprGroupPer10**
* RColorBrewer package of R: required by **agpdraw-line**

The wigToBigWig and bigWigToWig are available at 
http://hgdownload.cse.ucsc.edu/admin/exe/

You will also need *chrom.size* file to run **bam2bwshifted**. This file is a two-column text file of chromosome sizes. Please generate the file if you have not yet. In the case of human genome (hg19), simply run the following command after installing fetchChromSizes (UCSC):

```
fetchChromSizes hg19 > hg19.chrom.sizes
```

Usage
------------------

### agplus

generates a text table for the rapid and flexible creation of aggregation plots.

```
agplus [options] -b reference.bed [WIG/BEDGRAPH file(it expects STDIN if ommited)]
    -b, --reference=file             reference BED, requires strand(6th column) field (default: none, required)
    -d, --distance-from=origin       distance from start/center/end of the references is used (default: start)
    -o, --out=file                   output file name (default: STDOUT)
    -r, --range=from,to              aggregating range in base-pairs (default:-5000,5000)
    -a, --assignment=file            name->group assignment table, 2-column and tab-separated
```

### bam2bwshifted

converts BAM files into a coverage track (bigWig) of N bp shifted positions of the reads (fragment-midpoint). The midpoint counts are normalized as RPM (#reads per million reads).

```
bam2bwshifted [-o outfile.bw] [-s shiftsize (default: 73)] -g chrom.sizes BAM
```

### assignExprGroupPer10

generates *assignment* file of expression level groups (10%ile interval) from the *genes.fpkm_tracking* file of cufflinks (cuffdiff).

```
assignExprGroupPer10 [-o outputdir] [-c column-name-of-gene] genes.fpkm_tracking
```

### agpdraw-line

draws smoothed (Gaussian filtered) curves from the agplus output file. -s *sdev* sets the smoothing parameter of the Gaussian kernel (the bigger the *sdev*, the smoother the curve).

```
agpdraw-line [-o outfile.pdf] [-s sdev] [-r from,to] [-t title] agplus-output.txt
```

File formats
-------

### Target

The *target* of your ChIP-Seq analysis, should be in wiggle or bedgraph format.

Example:

```
chr1    9947    9948    0.052102
chr1    9952    9953    0.052102
chr1    9978    9979    0.052102
chr1    10003   10004   0.052102
chr1    10015   10016   0.052102
(...)
```

### Reference

The *reference* where you will aggregate *target* should be in BED6 format (6 columns BED). Typically, this file can be taken from a public database such as UCSC's table browser.

An example reference file of a gene locus looks like the following:

```
chr1    11873   14409   DDX11L1 0       +
chr1    14361   29370   WASH7P  0       -
chr1    34610   36081   FAM138F 0       -
chr1    34610   36081   FAM138A 0       -
chr1    69090   70008   OR4F5   0       +
(...)
```

### Assignment

The group *assignment* should be a two-column, tab-delimited text file. Typically, you might write the *assignment* definitions such as "gene-name[tab]group-name" per line in this file. Note that you do not need to write out all the names defined in *reference* here.

An example of *assignment* file of Up/Down/Stay genes groups looks like the following:

```
WASH7P	Upregulated
OR4F6	Upregulated
DDX11L1	Downregulated
FAM138A	Downregulated
FAM138F	Stay
(...)
```

### Output

The output file of **agplus** is simple tab-delimited table of average signal intensities in 1 bp per line. The table can be easily handled by subsequent analysis, using either gnuplot, R, Matlab or MS Excel etc.

Example:

```
distance	whole
-5000	0.0006606468404588093
-4999	0.0006476077580813323
-4998	0.0006345686757038559
(...)
4997    0.000534602377476536
4998    0.0004911387695516146
4999    0.0005128705735140753
5000    0.000554161001042751
```

Tutorial
--------

This tutorial aims to visualize H3K27ac ChIP-Seq signal distributions at transcription start sites (TSSs) grouped by the gene expression levels. The example data are H3K27ac ChIP-Seq of HeLa-S3 cells provided by the ENCODE/Broad Institute ([GSE29611](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29611)). The expression levels of HeLa-S3 were determined by ENCODE/Caltech RNA-Seq data ([GSE33480](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33480)).

### Prerequistics

You need to download the BAM file of H3K27ac ChIP-Seq from UCSC. 

* [ENCODE/Broad HeLa-S3 H3K27ac ChIP-Seq](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneHelas3H3k27acStdAlnRep1.bam)

We recommend using [UDR](http://rabadan.c2b2.columbia.edu/ENCODE/newsarch.html#091213) for faster downloading of ENCODE data from UCSC.

To define expression groups using **assignExprGroupPer10**, you need the *genes.fpkm_tracking* file of cufflinks (or cuffdiff) output using [this bam](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCaltechRnaSeq/wgEncodeCaltechRnaSeqHelas3R2x75Il200AlignsRep1V2.bam). For convenience, we prepared a pre-calculated *genes.fpkm_tracking* file of HeLa-S3. Please download [the file](http://chromatin.med.kyushu-u.ac.jp/maehara/agplus/genes.fpkm_tracking.gz) and decompress it.

Additionally, you need a gene locus definition file in BED6 format as *reference*. We prepared a parsed file from UCSC's hg19 refFlat definition for **agplus** (available from the table browser). Please download [the file](http://chromatin.med.kyushu-u.ac.jp/maehara/agplus/refFlat_hg19_simple.bed.gz) and decompress it.

### Four steps to generating an aggregation plot

Step 1: generate a coverage track (bigWig) file of fragment midpoints (`-s 100`; half-size of 200 bp) from the *target* BAM file

```
bam2bwshifted -s 100 -g hg19.chrom.sizes -o HeLaS3_H3K27ac.bw wgEncodeBroadHistoneHelas3H3k27acStdAlnRep1.bam
```

Step 2: decompress the bigWig file to treat with **agplus**

```
bigWigToWig HeLaS3_H3K27ac.bw HeLaS3_H3K27ac.wig
```

Step 3: aggregate *target* signals at all TSSs defined in *reference*

```
agplus -b refFlat_hg19_simple.bed -d start -o aggr_HeLaS3_H3K27ac.txt HeLaS3_H3K27ac.wig
```

The aggregation process takes about 2 minutes assuming the use of  Intel(R) Xeon(R) CPU E5-2687W @ 3.10GHz on x86_64 GNU/Linux)

>89.98s user 0.29s system 63% cpu 2:23.18 total

One can also combine Steps 2 & 3 to reduce the number of intermediate files.

```
bigWigToWig HeLaS3_H3K27ac.bw stdout | agplus -b refFlat_hg19_simple.bed -d start -o aggr_HeLaS3_H3K27ac.txt
```

Step4: draw the aggregation plot

```
agpdraw-line -t "H3K27ac at all TSSs" -r -2000,2000 -o H3K27acAllTSS.pdf aggr_HeLaS3_H3K27ac.txt
```

The result of a typical PDF file is shown below.

![H3K27ac at All TSSs](example/H3K27acAllTSS.png)

Please note that the aggregation is done against all TSSs when you omit the `-a` assignment file option, i.e. all gene TSSs are assigned to "whole" group.

### Grouping by assignment file

If you would like to see the distribution by expression levels, incorporate the following procedure after Step 2.

Create an assignment file from the genes.fpkm_tracking file.

```
assignExprGroupPer10 genes.fpkm_tracking
```

The output file is produced under the *assignedPer10* directory. Please check the file *egroup_ePer10.txt* that is automatically generated as the *assignment* file. The *assignment* file is the group definition of genes per 10%ile interval gene expression levels. (The output will be a multiple text file of different samples when the file is from *cuffdiff*.)

Then, Run agplus with the assignment file. 

```
bigWigToWig HeLaS3_H3K27ac.bw stdout | agplus -b refFlat_hg19_simple.bed -d start -a assignedPer10/egroup_ePer10.txt -o aggr_HeLaS3_H3K27ac_Expr10pt.txt
```

> 93.64s user 0.32s system 65% cpu 2:24.28 total

Finally, run agpdraw-line to get the PDF output.

```
agpdraw-line -t "H3K27ac by expression groups" -r -2000,2000 -o H3K27acByExpression.pdf aggr_HeLaS3_H3K27ac_Expr10pt.txt
```

The result:

![H3K27ac by Expression groups](example/H3K27acByExpression.png)


