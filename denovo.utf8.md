---
title: "Use of tximeta with a denovo transcriptome"
author: "Michael Love"
---

This is a project in collaboration with Lisa Johnson and Tessa Pierce
from the lab of Titus Brown.

# Transcript abundance

Transcript abundance was estimated with Salmon and the quantification
directories were collected as a single `.tar.gz` file and deposited to
Zenodo:

* https://zenodo.org/record/1486283#.W-xVBHpKiL8

These were downloaded and extracted to a directory
`Fundulus_rathbuni_salmon`, which contains 9 Salmon output
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

