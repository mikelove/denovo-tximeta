---
title: "Use of tximeta with a denovo transcriptome"
author: "Michael Love"
---

This is a project in collaboration with Lisa Johnson and Tessa Pierce
from the lab of Titus Brown.

# Transcript abundance

Transcript abundance was estimated with *Salmon* and the quantification
directories were collected as a single `.tar.gz` file and deposited to
Zenodo. It is important to upload the entire *Salmon* directories and
not just the `quant.sf` file, or else `tximeta` will not work.

* https://zenodo.org/record/1486283#.W-xVBHpKiL8

This file was downloaded and extracted to a directory
`Fundulus_rathbuni_salmon`, which contains 9 *Salmon* output
directories. 

# Attempting import without linking to annotation


```r
# install_github("mikelove/tximeta", dependencies=FALSE)
library(tximeta)
dir <- "Fundulus_rathbuni_salmon"
files <- file.path(list.files(dir, full=TRUE), "quant.sf")
all(file.exists(files))
```

```
## [1] TRUE
```

```r
names <- list.files(dir)
names <- sub(".quant","",names)
coldata <- data.frame(files, names)
se <- tximeta(coldata)
```

```
## importing quantifications
```

```
## reading in files with read_tsv
```

```
## 1 2 3 4 5 6 7 8 9 
## couldn't find matching transcriptome, returning un-ranged SummarizedExperiment
```

In the above code chunk, `tximeta` attempts to find a matching
transcriptome in its database but fails, and then returns an
"un-ranged SummarizedExperiment", meaning that we don't have any extra
information about the transcripts. We have the rownames and multiple
matrices of data:


```r
suppressPackageStartupMessages(library(SummarizedExperiment))
assayNames(se)
```

```
## [1] "counts"    "abundance" "length"
```

```r
assays(se)[["counts"]][1:4,1:4]
```

```
##                           F_rathbuni_BW_1 F_rathbuni_BW_2 F_rathbuni_BW_3
## TRINITY_DN114796_c0_g1_i1               4               0               4
## TRINITY_DN114738_c0_g1_i1               0               1               0
## TRINITY_DN114738_c1_g1_i1               1               0               0
## TRINITY_DN114736_c0_g1_i1               0               0               1
##                           F_rathbuni_FW_1
## TRINITY_DN114796_c0_g1_i1               6
## TRINITY_DN114738_c0_g1_i1               0
## TRINITY_DN114738_c1_g1_i1               1
## TRINITY_DN114736_c0_g1_i1               0
```

# Linking the quantifications to online annotations

Here, we have a de novo constructed transcriptome, and so if we want
to make the quantifications reproducible, we can formally link them to
their sequence (FASTA file) and additional information (GFF3 file),
which have also been deposited to Zenodo.

We provide the path to the *Salmon* index for the de novo
transcriptome, and provide links to Zenodo for the `fasta` and `gtf`
arguments. The `jsonFile` argument is the filename for the JSON file
that can be shared with others, linking these quantification files to
the appropriate metadata online.

We indicate that the source of the GFF3 file is the 
[dammit](http://dib-lab.github.io/dammit/) de novo transcriptome
annotator. This lets `tximeta` know that the *seqid* in the GFF3 file
should be used as the names of the transcripts (often the *seqid*
indicates the chromosome ID, but here we do not have chromosome
information).


```r
makeLinkedTxome(indexDir="F_rathbuni.trinity_out",
                source="dammit",
                organism="Fundulus rathbuni",
                release="0",
                genome="none",
                fasta="https://zenodo.org/record/1486276/files/F_rathbuni.trinity_out.fasta",
                gtf="https://zenodo.org/record/2226742/files/F_rathbuni.trinity_out.Trinity.fasta.dammit.gff3",
                jsonFile="F_rathbuni.json")
```


```
## writing linkedTxome to F_rathbuni.json
```

```
## saving linkedTxome in bfc
```

Now we can try importing the quantifications again:


```r
se <- tximeta(coldata)
```

```
## importing quantifications
```

```
## reading in files with read_tsv
```

```
## 1 2 3 4 5 6 7 8 9 
## found matching linked transcriptome:
## [ dammit - Fundulus rathbuni - release 0 ]
## loading existing TxDb created: 2018-12-13 18:26:20
## Loading required package: GenomicFeatures
## Loading required package: AnnotationDbi
## generating transcript ranges
```

```
## Warning in checkAssays2Txps(assays, txps): missing some transcripts!
##  357723 out of 501215 are missing from the GTF and dropped from SummarizedExperiment output
```

This time, `tximeta` was able to recognize the transcriptome, and to
pull down the appropriate metadata from Zenodo. The quantifications
are also now linked to the appropriate transcript FASTA file, enabling
computational reproducibility of the quantification step.

The warning about missing ranges from the annotation file is because
the Bioconductor GTF/GFF importer uses `mRNA` as an indicator for
transcripts, and not all of the transcripts in the quantification
files have corresponding `mRNA` entries in the GFF3 file. This is an
area for development...

We can see the additional information attached to the rows of our
quantification data:


```r
rowRanges(se)
```

```
## GRanges object with 143492 ranges and 3 metadata columns:
##                                              seqnames    ranges strand |
##                                                 <Rle> <IRanges>  <Rle> |
##   TRINITY_DN114791_c0_g1_i1 TRINITY_DN114791_c0_g1_i1    1-2308      + |
##   TRINITY_DN114724_c0_g2_i1 TRINITY_DN114724_c0_g2_i1     1-635      - |
##   TRINITY_DN114714_c0_g1_i1 TRINITY_DN114714_c0_g1_i1     1-282      - |
##   TRINITY_DN114702_c1_g1_i1 TRINITY_DN114702_c1_g1_i1     1-818      + |
##   TRINITY_DN114789_c0_g1_i1 TRINITY_DN114789_c0_g1_i1     1-362      - |
##                         ...                       ...       ...    ... .
##    TRINITY_DN82378_c0_g1_i1  TRINITY_DN82378_c0_g1_i1     1-698      - |
##    TRINITY_DN82373_c0_g1_i1  TRINITY_DN82373_c0_g1_i1     1-259      + |
##    TRINITY_DN89144_c0_g1_i1  TRINITY_DN89144_c0_g1_i1     1-558      - |
##    TRINITY_DN70896_c0_g1_i1  TRINITY_DN70896_c0_g1_i1     1-322      + |
##    TRINITY_DN70859_c1_g1_i1  TRINITY_DN70859_c1_g1_i1     1-557      - |
##                                 tx_id
##                             <integer>
##   TRINITY_DN114791_c0_g1_i1      1290
##   TRINITY_DN114724_c0_g2_i1      1283
##   TRINITY_DN114714_c0_g1_i1      1282
##   TRINITY_DN114702_c1_g1_i1      1280
##   TRINITY_DN114789_c0_g1_i1      1289
##                         ...       ...
##    TRINITY_DN82378_c0_g1_i1    204616
##    TRINITY_DN82373_c0_g1_i1    204615
##    TRINITY_DN89144_c0_g1_i1    205073
##    TRINITY_DN70896_c0_g1_i1    204113
##    TRINITY_DN70859_c1_g1_i1    204112
##                                                                                                               gene_id
##                                                                                                       <CharacterList>
##   TRINITY_DN114791_c0_g1_i1                           ORF Transcript_5|g.2 Transcript_5|m.2 type:complete len:190 (+)
##   TRINITY_DN114724_c0_g2_i1                    ORF Transcript_13|g.8 Transcript_13|m.8 type:5prime_partial len:83 (-)
##   TRINITY_DN114714_c0_g1_i1                        ORF Transcript_36|g.13 Transcript_36|m.13 type:internal len:95 (-)
##   TRINITY_DN114702_c1_g1_i1                        ORF Transcript_43|g.18 Transcript_43|m.18 type:complete len:82 (+)
##   TRINITY_DN114789_c0_g1_i1                       ORF Transcript_44|g.19 Transcript_44|m.19 type:internal len:121 (-)
##                         ...                                                                                       ...
##    TRINITY_DN82378_c0_g1_i1        ORF Transcript_500975|g.496894 Transcript_500975|m.496894 type:complete len:94 (-)
##    TRINITY_DN82373_c0_g1_i1        ORF Transcript_500979|g.496896 Transcript_500979|m.496896 type:internal len:87 (+)
##    TRINITY_DN89144_c0_g1_i1  ORF Transcript_501115|g.496901 Transcript_501115|m.496901 type:3prime_partial len:97 (-)
##    TRINITY_DN70896_c0_g1_i1 ORF Transcript_501179|g.496904 Transcript_501179|m.496904 type:5prime_partial len:105 (+)
##    TRINITY_DN70859_c1_g1_i1  ORF Transcript_501212|g.496910 Transcript_501212|m.496910 type:5prime_partial len:94 (-)
##                                                                                                               tx_name
##                                                                                                           <character>
##   TRINITY_DN114791_c0_g1_i1                           ORF Transcript_5|g.2 Transcript_5|m.2 type:complete len:190 (+)
##   TRINITY_DN114724_c0_g2_i1                    ORF Transcript_13|g.8 Transcript_13|m.8 type:5prime_partial len:83 (-)
##   TRINITY_DN114714_c0_g1_i1                        ORF Transcript_36|g.13 Transcript_36|m.13 type:internal len:95 (-)
##   TRINITY_DN114702_c1_g1_i1                        ORF Transcript_43|g.18 Transcript_43|m.18 type:complete len:82 (+)
##   TRINITY_DN114789_c0_g1_i1                       ORF Transcript_44|g.19 Transcript_44|m.19 type:internal len:121 (-)
##                         ...                                                                                       ...
##    TRINITY_DN82378_c0_g1_i1        ORF Transcript_500975|g.496894 Transcript_500975|m.496894 type:complete len:94 (-)
##    TRINITY_DN82373_c0_g1_i1        ORF Transcript_500979|g.496896 Transcript_500979|m.496896 type:internal len:87 (+)
##    TRINITY_DN89144_c0_g1_i1  ORF Transcript_501115|g.496901 Transcript_501115|m.496901 type:3prime_partial len:97 (-)
##    TRINITY_DN70896_c0_g1_i1 ORF Transcript_501179|g.496904 Transcript_501179|m.496904 type:5prime_partial len:105 (+)
##    TRINITY_DN70859_c1_g1_i1  ORF Transcript_501212|g.496910 Transcript_501212|m.496910 type:5prime_partial len:94 (-)
##   -------
##   seqinfo: 143492 sequences from an unspecified genome; no seqlengths
```




```r
session_info()
```

```
## Session info -------------------------------------------------------------
```

```
##  setting  value                                             
##  version  R Under development (unstable) (2018-05-14 r74725)
##  system   x86_64, darwin15.6.0                              
##  ui       X11                                               
##  language (EN)                                              
##  collate  en_US.UTF-8                                       
##  tz       America/New_York                                  
##  date     2018-12-13
```

```
## Packages -----------------------------------------------------------------
```

```
##  package              * version   date       source        
##  AnnotationDbi        * 1.43.1    2018-06-13 Bioconductor  
##  AnnotationFilter       1.5.2     2018-06-13 Bioconductor  
##  assertthat             0.2.0     2017-04-11 CRAN (R 3.6.0)
##  backports              1.1.2     2017-12-13 CRAN (R 3.6.0)
##  base                 * 3.6.0     2018-05-15 local         
##  bindr                  0.1.1     2018-03-13 CRAN (R 3.6.0)
##  bindrcpp             * 0.2.2     2018-03-29 CRAN (R 3.6.0)
##  Biobase              * 2.41.1    2018-06-17 Bioconductor  
##  BiocFileCache        * 1.5.5     2018-07-20 Bioconductor  
##  BiocGenerics         * 0.27.1    2018-06-17 Bioconductor  
##  BiocParallel         * 1.15.6    2018-06-28 Bioconductor  
##  biomaRt                2.37.3    2018-06-29 Bioconductor  
##  Biostrings             2.49.0    2018-05-22 Bioconductor  
##  bit                    1.1-14    2018-05-29 CRAN (R 3.6.0)
##  bit64                  0.9-7     2017-05-08 CRAN (R 3.6.0)
##  bitops                 1.0-6     2013-08-17 CRAN (R 3.6.0)
##  blob                   1.1.1     2018-03-25 CRAN (R 3.6.0)
##  compiler               3.6.0     2018-05-15 local         
##  crayon                 1.3.4     2017-09-16 CRAN (R 3.6.0)
##  curl                   3.2       2018-03-28 CRAN (R 3.6.0)
##  datasets             * 3.6.0     2018-05-15 local         
##  DBI                    1.0.0     2018-05-02 CRAN (R 3.6.0)
##  dbplyr               * 1.2.1     2018-02-19 CRAN (R 3.6.0)
##  DelayedArray         * 0.7.30    2018-08-17 Bioconductor  
##  devtools             * 1.13.6    2018-06-27 CRAN (R 3.6.0)
##  digest                 0.6.15    2018-01-28 CRAN (R 3.6.0)
##  dplyr                  0.7.6     2018-06-29 CRAN (R 3.6.0)
##  ensembldb              2.5.8     2018-08-31 Bioconductor  
##  evaluate               0.10.1    2017-06-24 CRAN (R 3.6.0)
##  GenomeInfoDb         * 1.17.1    2018-06-13 Bioconductor  
##  GenomeInfoDbData       1.1.0     2018-05-15 Bioconductor  
##  GenomicAlignments      1.17.2    2018-06-13 Bioconductor  
##  GenomicFeatures      * 1.33.2    2018-08-13 Bioconductor  
##  GenomicRanges        * 1.33.6    2018-06-13 Bioconductor  
##  glue                   1.2.0     2017-10-29 CRAN (R 3.6.0)
##  graphics             * 3.6.0     2018-05-15 local         
##  grDevices            * 3.6.0     2018-05-15 local         
##  grid                   3.6.0     2018-05-15 local         
##  hms                    0.4.2     2018-03-10 CRAN (R 3.6.0)
##  htmltools              0.3.6     2017-04-28 CRAN (R 3.6.0)
##  httr                   1.3.1     2017-08-20 CRAN (R 3.6.0)
##  IRanges              * 2.15.14   2018-06-13 Bioconductor  
##  jsonlite               1.5       2017-06-01 CRAN (R 3.6.0)
##  knitr                  1.20      2018-02-20 CRAN (R 3.6.0)
##  lattice                0.20-35   2017-03-25 CRAN (R 3.6.0)
##  lazyeval               0.2.1     2017-10-29 CRAN (R 3.6.0)
##  magrittr               1.5       2014-11-22 CRAN (R 3.6.0)
##  Matrix                 1.2-14    2018-04-13 CRAN (R 3.6.0)
##  matrixStats          * 0.54.0    2018-07-23 cran (@0.54.0)
##  memoise                1.1.0     2017-04-21 CRAN (R 3.6.0)
##  methods              * 3.6.0     2018-05-15 local         
##  parallel             * 3.6.0     2018-05-15 local         
##  pillar                 1.2.3     2018-05-25 CRAN (R 3.6.0)
##  pkgconfig              2.0.1     2017-03-21 CRAN (R 3.6.0)
##  prettyunits            1.0.2     2015-07-13 CRAN (R 3.6.0)
##  progress               1.2.0     2018-06-14 CRAN (R 3.6.0)
##  ProtGenerics           1.13.0    2018-06-13 Bioconductor  
##  purrr                  0.2.5     2018-05-29 CRAN (R 3.6.0)
##  R6                     2.2.2     2017-06-17 CRAN (R 3.6.0)
##  rappdirs               0.3.1     2016-03-28 CRAN (R 3.6.0)
##  Rcpp                   0.12.19   2018-10-01 CRAN (R 3.6.0)
##  RCurl                  1.95-4.10 2018-01-04 CRAN (R 3.6.0)
##  readr                  1.1.1     2017-05-16 CRAN (R 3.6.0)
##  rlang                  0.2.0     2018-02-20 CRAN (R 3.6.0)
##  rmarkdown            * 1.10      2018-06-11 CRAN (R 3.6.0)
##  rprojroot              1.3-2     2018-01-03 CRAN (R 3.6.0)
##  Rsamtools              1.33.7    2018-10-19 Bioconductor  
##  RSQLite                2.1.1     2018-05-06 CRAN (R 3.6.0)
##  rtracklayer            1.41.3    2018-06-13 Bioconductor  
##  S4Vectors            * 0.19.19   2018-07-18 Bioconductor  
##  stats                * 3.6.0     2018-05-15 local         
##  stats4               * 3.6.0     2018-05-15 local         
##  stringi                1.2.2     2018-05-02 CRAN (R 3.6.0)
##  stringr                1.3.1     2018-05-10 CRAN (R 3.6.0)
##  SummarizedExperiment * 1.11.6    2018-07-17 Bioconductor  
##  testthat             * 2.0.0     2017-12-13 CRAN (R 3.6.0)
##  tibble                 1.4.2     2018-01-22 CRAN (R 3.6.0)
##  tidyselect             0.2.4     2018-02-26 CRAN (R 3.6.0)
##  tools                  3.6.0     2018-05-15 local         
##  tximeta              * 1.1.9     2018-12-13 Bioconductor  
##  tximport               1.11.4    2018-11-27 Bioconductor  
##  utils                * 3.6.0     2018-05-15 local         
##  withr                  2.1.2     2018-03-15 CRAN (R 3.6.0)
##  XML                    3.98-1.11 2018-04-16 CRAN (R 3.6.0)
##  XVector                0.21.3    2018-06-23 Bioconductor  
##  yaml                   2.1.19    2018-05-01 CRAN (R 3.6.0)
##  zlibbioc               1.27.0    2018-05-15 Bioconductor
```
