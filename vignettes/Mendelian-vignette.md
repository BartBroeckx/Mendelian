---
title: "Mendelian vignette"
author: "Bart Broeckx"
date: "2015-10-27"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Mendelian-vignette}
  %\VignetteEngine{knitr::rmarkdown}
  %\usepackage[utf8]{inputenc}
---
# Introduction

Since the development of second generation sequencing, the potential to sequence large numbers of individuals has tremendously increased, with, at the same time, a decrease in costs. To analyze the vast amounts of data generated, a wide variety of algorithms and pipelines have been developped for e.g. read mapping and variant calling. However, in disease association studies, the real search only starts after the variant calling. The identification of variant(s) or gene(s) responsible for a (Mendelian or complex) disorder, can pose quite a challenge as variant calling often results in a large amount of variants and ideally, we should be able to reduce this to a very limited number of variants or genes (1?). This package aims to be the link between variant calling and assessing the effect of mutations on protein functionality (by e.g. [PolyPhen](http://genetics.bwh.harvard.edu/pph2/) or [Provean](http://provean.jcvi.org/index.php)).   

The Mendelian package was designed to identify causal variants for Mendelian disorders, assuming a dominant or recessive mode of inheritance and allows for variable degrees of penetrance, detectance and allelic heterogeneity. It can be used to process both the standard Variant Call Format (VCF) and the output from the [CLC genomics workbench] (http://www.clcbio.com/products/clc-genomics-workbench/). Filtering variants present in popular database (e.g. dbSNP) is also supported.  

This vignette is divided in three large parts:
* Preprocessing
* Variant filtering
* Combining filters

The sequence of functions can also be found in the following graph:
![](http://www.heupdysplasie.ugent.be/Mendelian/figure2.png)

*varfilter* is not mentioned as it can be used anytime, with one exception: VCF files should first be processed with *VCFfile* before *varfilter* is used.

# Preprocessing
This package supports two commonly used formats containing sequencing variants:
* VCF
* txt output from the CLC genomics workbench

VCF files consist of a header and a variant section. The variant section itself has 8 mandatory columns. More information on VCF files can be found [here] (http://bioinformatics.oxfordjournals.org/content/27/15/2156). Output from the CLC genomics workbench consists of a variable amount of columns, depending on the processing steps done inside CLC. As the aim is reducing the number of variants, I personally recommend to export the files after only nonsynonimous variants are retained.

Depending on the filters you wish to use downstream, the CLC files should be preprocessed or not. VCF files however, should always be preprocessed by the function *VCFfile*. 

## VCFfile
*VCFfile* allows the one by one preprocessing of VCF files. The main goal is the selection of the non-reference variants for a certain individual and reorganising the file to make it compatible with the filter functions. In addition, a filter can be specified, retaining only those variants with e.g. a certain quality score. As VCF files might contain variant information of several individuals, the column name containing variant information of the individual of interest should always be specified. 

The "test" data frame has 21 rows, with row 12 (chr1 1581713) containing 2 non-reference variants. 

```
## Loading required package: stringr
```


```r
library("Mendelian")
```

```
## 
## Attaching package: 'Mendelian'
## 
## The following object is masked _by_ '.GlobalEnv':
## 
##     VCFfile
```

```r
  dim(test)
```

```
## [1] 21 10
```

```r
  str(test)
```

```
## 'data.frame':	21 obs. of  10 variables:
##  $ V1 : Factor w/ 24 levels "chr1","chr10",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ V2 : int  1563789 1564952 1570039 1575616 1577003 1580738 1581065 1581096 1581558 1581586 ...
##  $ V3 : Factor w/ 54924 levels ".","rs10000250",..: 47593 4922 1 51334 1 1 1 1 8702 40635 ...
##  $ V4 : Factor w/ 693 levels "A","AAAAAAAAAAAAAAAAAAG",..: 176 609 367 527 176 527 367 176 176 527 ...
##  $ V5 : Factor w/ 634 levels "A","A,AC","A,ACTGCTGCTG",..: 1 468 152 152 468 1 152 468 468 152 ...
##  $ V6 : num  946 529 101 1530 119 ...
##  $ V7 : Factor w/ 80 levels "HARD_TO_VALIDATE",..: 69 80 31 30 31 80 79 79 31 31 ...
##  $ V8 : Factor w/ 61487 levels "AB=0.004;AC=1;AF=0.50;AN=2;BaseQRankSum=1.605;DB;DP=480;Dels=0.00;FS=8.451;HRun=0;HaplotypeScore=2.6000;MQ=4.51;MQ0=449;MQRankS"| __truncated__,..: 41114 45274 60986 58527 61394 28728 17694 16338 54678 53041 ...
##  $ V9 : Factor w/ 1 level "GT:AD:DP:GQ:PL": 1 1 1 1 1 1 1 1 1 1 ...
##  $ V10: Factor w/ 60111 levels "0/1:1,10:11:6.90:348,0,7",..: 55487 55182 54866 52753 58960 27959 17315 22508 49617 48313 ...
```

```r
  test[,c(1,2,4,5,9)]
```

```
##       V1      V2 V4  V5             V9
## 160 chr1 1563789  C   A GT:AD:DP:GQ:PL
## 161 chr1 1564952 TG   T GT:AD:DP:GQ:PL
## 162 chr1 1570039  G   C GT:AD:DP:GQ:PL
## 163 chr1 1575616  T   C GT:AD:DP:GQ:PL
## 164 chr1 1577003  C   T GT:AD:DP:GQ:PL
## 165 chr1 1580738  T   A GT:AD:DP:GQ:PL
## 166 chr1 1581065  G   C GT:AD:DP:GQ:PL
## 167 chr1 1581096  C   T GT:AD:DP:GQ:PL
## 168 chr1 1581558  C   T GT:AD:DP:GQ:PL
## 169 chr1 1581586  T   C GT:AD:DP:GQ:PL
## 170 chr1 1581604  G   A GT:AD:DP:GQ:PL
## 171 chr1 1581713  A C,G GT:AD:DP:GQ:PL
## 172 chr1 1581759  G   A GT:AD:DP:GQ:PL
## 173 chr1 1581881  T   C GT:AD:DP:GQ:PL
## 174 chr1 1582106  T   C GT:AD:DP:GQ:PL
## 175 chr1 1582202  C   T GT:AD:DP:GQ:PL
## 176 chr1 1585597  A   G GT:AD:DP:GQ:PL
## 177 chr1 1585642  G   T GT:AD:DP:GQ:PL
## 178 chr1 1586752  T   C GT:AD:DP:GQ:PL
## 179 chr1 1588670  G   A GT:AD:DP:GQ:PL
## 180 chr1 1588717  G   A GT:AD:DP:GQ:PL
```


With the filter, the number of variants is reduced to 6, with one region (1581713) containing two non-reference variants. After preprocessing without a filter, this
renders all the variants, but with 22 rows now as region 1581713 contains two non-reference variants.

```r
# test 1: filter = TRUE
  x <- test
  sample <- "V10"
  filter <- TRUE
  value <- "PASS"
  VCFfile(x,sample, filter, value)
```

```
##      Chromosome  Region Allele     Zygosity
## 160        chr1 1563789      A   Homozygous
## 171        chr1 1581713      C Heterozygous
## 1711       chr1 1581713      G Heterozygous
## 176        chr1 1585597      G   Homozygous
## 177        chr1 1585642      T   Homozygous
## 178        chr1 1586752      C   Homozygous
```

```r
# test 2: filter= FALSE
  x <- test
  sample <- "V10"
  filter <- FALSE
  VCFfile(x,sample, filter)
```

```
##      Chromosome  Region Allele     Zygosity
## 160        chr1 1563789      A   Homozygous
## 161        chr1 1564952      T   Homozygous
## 162        chr1 1570039      C   Homozygous
## 163        chr1 1575616      C   Homozygous
## 164        chr1 1577003      T   Homozygous
## 165        chr1 1580738      A Heterozygous
## 166        chr1 1581065      C Heterozygous
## 167        chr1 1581096      T Heterozygous
## 168        chr1 1581558      T   Homozygous
## 169        chr1 1581586      C   Homozygous
## 170        chr1 1581604      A Heterozygous
## 171        chr1 1581713      C Heterozygous
## 1711       chr1 1581713      G Heterozygous
## 172        chr1 1581759      A Heterozygous
## 173        chr1 1581881      C Heterozygous
## 174        chr1 1582106      C   Homozygous
## 175        chr1 1582202      T Heterozygous
## 176        chr1 1585597      G   Homozygous
## 177        chr1 1585642      T   Homozygous
## 178        chr1 1586752      C   Homozygous
## 179        chr1 1588670      A Heterozygous
## 180        chr1 1588717      A Heterozygous
```

## CLCfile
*CLCfile* allows the one by one preprocessing of output from the CLC genomics workbench. I personally recommend to export only after filtering for nonsynonymous variants has been done, as this 1) reduces the number of variants 2) adds a functional unit (e.g. a gene) that can be used for filtering with gDom and gRec. 

!! IMPORTANT: the function CLCfile should NOT be used prior to nDom and nRec !!


```r
str(CLCfile1)
```

```
## 'data.frame':	8 obs. of  21 variables:
##  $ Chromosome                 : chr  "chr1" "chr1" "chr1" "chr1" ...
##  $ Region                     : chr  "266523" "363618" "604894" "699873" ...
##  $ Type                       : chr  "SNV" "SNV" "SNV" "SNV" ...
##  $ Reference                  : chr  "G" "A" "A" "G" ...
##  $ Allele                     : chr  "T" "G" "G" "A" ...
##  $ Reference.allele           : chr  "No" "No" "No" "No" ...
##  $ Length                     : int  1 1 1 1 1 1 1 3
##  $ Linkage                    : logi  NA NA NA NA NA NA ...
##  $ Zygosity                   : chr  "Heterozygous" "Heterozygous" "Homozygous" "Homozygous" ...
##  $ Count                      : int  12 27 212 21 23 15 22 5
##  $ Coverage                   : int  24 52 212 21 23 61 62 5
##  $ Frequency                  : num  50 51.9 100 100 100 ...
##  $ Probability                : num  1 1 1 1 1 1 1 1
##  $ Forward.read.count         : int  8 27 134 17 23 7 5 1
##  $ Reverse.read.count         : int  5 3 87 7 2 9 17 4
##  $ Forward.reverse.balance    : num  0.385 0.1 0.394 0.292 0.08 0.438 0.227 0.2
##  $ Average.quality            : num  38 36.1 36.2 35.1 33 ...
##  $ Coding.region.change       : chr  "ENSCAFT00000000001:c.1802C>A" "ENSCAFT00000000003:c.12A>G" "ENSCAFT00000000008:c.775T>C" "ENSCAFT00000037436:c.371C>T" ...
##  $ Amino.acid.change          : chr  "ENSCAFT00000000001:p.Pro601His" "ENSCAFT00000000003:p.Ile4Met" "ENSCAFT00000000008:p.Phe259Leu" "ENSCAFT00000037436:p.Ala124Val" ...
##  $ Other.variants.within.codon: chr  "No" "No" "No" "No" ...
##  $ Non.synonymous             : chr  "Yes" "Yes" "Yes" "Yes" ...
```
All three objects contain 21 columns and 8 or 9 rows. Each row has columns specifying the chromosome, position ("region") and associated variant ("allele") found. 


As we might expect, one variant might be associated with more than one functional unit. This is demonstrated at line 7 in CLCfile1 where the variant causes an amino-acid change in both ENSCAFT00000043665 and ENSCAFT00000048679:


```r
CLCfile1[7,19]
```

```
## [1] "ENSCAFT00000043665:p.[Arg194Gln]; ENSCAFT00000048679:p.[Ser57Leu]"
```
If preprocessing is not done with the CLCfile function, this will be an issue later on in gDom and gRec as only the first functional unit will be kept for filtering. The difference between specifying multiple as TRUE and FALSE is demonstrated in the following two examples: 


```r
diff <- CLCfile1[7,]
# example 1
CLCfile1proc <-CLCfile(diff, "Coding.region.change", TRUE, "; ")
nrow(CLCfile1proc)
```

```
## [1] 2
```

```r
CLCfile1proc
```

```
##      Chromosome   Region Type Reference Allele Reference.allele Length
## 23         chr1 13397669  SNV         C      T               No      1
## 23.1       chr1 13397669  SNV         C      T               No      1
##      Linkage     Zygosity Count Coverage Frequency Probability
## 23        NA Heterozygous    22       62    35.484           1
## 23.1      NA Heterozygous    22       62    35.484           1
##      Forward.read.count Reverse.read.count Forward.reverse.balance
## 23                    5                 17                   0.227
## 23.1                  5                 17                   0.227
##      Average.quality
## 23            34.591
## 23.1          34.591
##                                              Coding.region.change
## 23   ENSCAFT00000043665:c.[581G>A]; ENSCAFT00000048679:c.[170C>T]
## 23.1 ENSCAFT00000043665:c.[581G>A]; ENSCAFT00000048679:c.[170C>T]
##                                                      Amino.acid.change
## 23   ENSCAFT00000043665:p.[Arg194Gln]; ENSCAFT00000048679:p.[Ser57Leu]
## 23.1 ENSCAFT00000043665:p.[Arg194Gln]; ENSCAFT00000048679:p.[Ser57Leu]
##      Other.variants.within.codon Non.synonymous              group
## 23                           Yes            Yes ENSCAFT00000043665
## 23.1                         Yes            Yes ENSCAFT00000048679
```

```r
# example 2
CLCfile1proc <-CLCfile(diff, "Coding.region.change", FALSE)
nrow(CLCfile1proc)
```

```
## [1] 1
```

```r
CLCfile1proc
```

```
##    Chromosome   Region Type Reference Allele Reference.allele Length
## 23       chr1 13397669  SNV         C      T               No      1
##    Linkage     Zygosity Count Coverage Frequency Probability
## 23      NA Heterozygous    22       62    35.484           1
##    Forward.read.count Reverse.read.count Forward.reverse.balance
## 23                  5                 17                   0.227
##    Average.quality
## 23          34.591
##                                            Coding.region.change
## 23 ENSCAFT00000043665:c.[581G>A]; ENSCAFT00000048679:c.[170C>T]
##                                                    Amino.acid.change
## 23 ENSCAFT00000043665:p.[Arg194Gln]; ENSCAFT00000048679:p.[Ser57Leu]
##    Other.variants.within.codon Non.synonymous              group
## 23                         Yes            Yes ENSCAFT00000043665
```
As clearly demonstrated, the first example results in 2 rows, each containing one functional unit. The second example only retains the first gene, ENSCAFT00000043665.

In the previous example, the column "Coding.region.change" was used as the reference column containing the functional unit. However, a different column ("Amino.acid.change") can also be used. Both columns are default output from the CLC Genomics Workbench if a check for nonsynonimous variants is done. Chosing a different column ("Amino.acid.change" instead of "Coding.region.change"), should give the same result after processing:

```r
diff <- CLCfile1[7,]
# example 1
CLCfile1proc <-CLCfile(diff, "Coding.region.change", TRUE, "; ")
nrow(CLCfile1proc)
```

```
## [1] 2
```

```r
CLCfile1proc[,c(1,18,19)]
```

```
##      Chromosome
## 23         chr1
## 23.1       chr1
##                                              Coding.region.change
## 23   ENSCAFT00000043665:c.[581G>A]; ENSCAFT00000048679:c.[170C>T]
## 23.1 ENSCAFT00000043665:c.[581G>A]; ENSCAFT00000048679:c.[170C>T]
##                                                      Amino.acid.change
## 23   ENSCAFT00000043665:p.[Arg194Gln]; ENSCAFT00000048679:p.[Ser57Leu]
## 23.1 ENSCAFT00000043665:p.[Arg194Gln]; ENSCAFT00000048679:p.[Ser57Leu]
```

```r
# example 2
CLCfile1proc <-CLCfile(diff, "Amino.acid.change", TRUE, "; ")
nrow(CLCfile1proc)
```

```
## [1] 2
```

```r
CLCfile1proc[,c(1,18,19)]
```

```
##      Chromosome
## 23         chr1
## 23.1       chr1
##                                              Coding.region.change
## 23   ENSCAFT00000043665:c.[581G>A]; ENSCAFT00000048679:c.[170C>T]
## 23.1 ENSCAFT00000043665:c.[581G>A]; ENSCAFT00000048679:c.[170C>T]
##                                                      Amino.acid.change
## 23   ENSCAFT00000043665:p.[Arg194Gln]; ENSCAFT00000048679:p.[Ser57Leu]
## 23.1 ENSCAFT00000043665:p.[Arg194Gln]; ENSCAFT00000048679:p.[Ser57Leu]
```

## annot
The *annot* function can be used to annotate VCF files preprocessed by *VCFfile*. 
Output from the CLC genomics workbench that is not annotated, should be processed immediately by this function instead of using the CLCfile function. Again, please note that both annot and CLCfile should not be used BEFORE nDom and nRec are used. 

Using *annot* on CLCfile1 yields the following result:




```r
# produces an annotated CLC file with all variants retained.
AnnotCLCfile1 <- annot(CLCfile1, genBED, type="BED", nomatch=TRUE, CLC=TRUE)
nrow(AnnotCLCfile1)
```

```
## [1] 10
```

```r
nrow(CLCfile1)
```

```
## [1] 8
```

```r
AnnotCLCfile1[,c(1,2,18,22)]
```

```
##      Chromosome             Region
## 1          chr1             266523
## 2          chr1             363618
## 3          chr1             604894
## 4          chr1             699873
## 5          chr1      758487^758488
## 22         chr1           13397637
## 22.1       chr1           13397637
## 23         chr1           13397669
## 23.1       chr1           13397669
## 48         chr1 16806440..16806442
##                                              Coding.region.change
## 1                                    ENSCAFT00000000001:c.1802C>A
## 2                                      ENSCAFT00000000003:c.12A>G
## 3                                     ENSCAFT00000000008:c.775T>C
## 4                                     ENSCAFT00000037436:c.371C>T
## 5                                  ENSCAFT00000000011:c.74_75insG
## 22   ENSCAFT00000043665:c.[613A>T]; ENSCAFT00000048679:c.[138T>A]
## 22.1 ENSCAFT00000043665:c.[613A>T]; ENSCAFT00000048679:c.[138T>A]
## 23   ENSCAFT00000043665:c.[581G>A]; ENSCAFT00000048679:c.[170C>T]
## 23.1 ENSCAFT00000043665:c.[581G>A]; ENSCAFT00000048679:c.[170C>T]
## 48                             ENSCAFT00000000149:c.565_567delGAG
##              annotation
## 1    ENSCAFT00000000001
## 2    ENSCAFT00000000003
## 3    ENSCAFT00000000008
## 4    ENSCAFT00000037436
## 5    ENSCAFT00000000011
## 22   ENSCAFT00000043665
## 22.1 ENSCAFT00000048679
## 23   ENSCAFT00000043665
## 23.1 ENSCAFT00000048679
## 48   ENSCAFT00000000149
```

```r
CLCfile1[,c(1,2,18)]
```

```
##    Chromosome             Region
## 1        chr1             266523
## 2        chr1             363618
## 3        chr1             604894
## 4        chr1             699873
## 5        chr1      758487^758488
## 22       chr1           13397637
## 23       chr1           13397669
## 48       chr1 16806440..16806442
##                                            Coding.region.change
## 1                                  ENSCAFT00000000001:c.1802C>A
## 2                                    ENSCAFT00000000003:c.12A>G
## 3                                   ENSCAFT00000000008:c.775T>C
## 4                                   ENSCAFT00000037436:c.371C>T
## 5                                ENSCAFT00000000011:c.74_75insG
## 22 ENSCAFT00000043665:c.[613A>T]; ENSCAFT00000048679:c.[138T>A]
## 23 ENSCAFT00000043665:c.[581G>A]; ENSCAFT00000048679:c.[170C>T]
## 48                           ENSCAFT00000000149:c.565_567delGAG
```

As demonstrated, the annotation column added by *annot* yields identical results as the columns "Coding.region.change" or "Amino.acid.change". However, this is not necessarily the case. *Annot* does not take aminoacid changes into account, while both other columns do. It might be that a variant lies within the borders of an additional gene, but has no effect on the aminoacid composition of the protein. This will result in an additional row when *annot* is used. When comparing "Coding.region.change" or "Amino.acid.change" with the result of annot, annot should always result in at least the same number of rows (but might give more).

In addition, although the previous example uses a "BED file", a "GTF" can also be used.

An additional option specifies what should be done with variants that do not lie within any functional units. Nomatch should be TRUE when all variants should be kept, it should be FALSE when only variants within the functional unit should be kept. To demonstrate the effect of both options, we change the region of the first variant. 


```r
CLCfile1[1,2]
```

```
## [1] "266523"
```

```r
a <- CLCfile1
a[1,2] <- "240000" 
a[c(1,2),]
```

```
##   Chromosome Region Type Reference Allele Reference.allele Length Linkage
## 1       chr1 240000  SNV         G      T               No      1      NA
## 2       chr1 363618  SNV         A      G               No      1      NA
##       Zygosity Count Coverage Frequency Probability Forward.read.count
## 1 Heterozygous    12       24    50.000           1                  8
## 2 Heterozygous    27       52    51.923           1                 27
##   Reverse.read.count Forward.reverse.balance Average.quality
## 1                  5                   0.385          38.000
## 2                  3                   0.100          36.067
##           Coding.region.change              Amino.acid.change
## 1 ENSCAFT00000000001:c.1802C>A ENSCAFT00000000001:p.Pro601His
## 2   ENSCAFT00000000003:c.12A>G   ENSCAFT00000000003:p.Ile4Met
##   Other.variants.within.codon Non.synonymous
## 1                          No            Yes
## 2                          No            Yes
```

```r
# produces an annotated CLC file with all variants retained.
 
AnnotCLCfile1 <- annot(a, genBED, type="BED", nomatch=TRUE, CLC=TRUE)
nrow(AnnotCLCfile1)
```

```
## [1] 10
```

```r
AnnotCLCfile1[c(1:2),]
```

```
##   Chromosome Region Type Reference Allele Reference.allele Length Linkage
## 1       chr1 240000  SNV         G      T               No      1      NA
## 2       chr1 363618  SNV         A      G               No      1      NA
##       Zygosity Count Coverage Frequency Probability Forward.read.count
## 1 Heterozygous    12       24    50.000           1                  8
## 2 Heterozygous    27       52    51.923           1                 27
##   Reverse.read.count Forward.reverse.balance Average.quality
## 1                  5                   0.385          38.000
## 2                  3                   0.100          36.067
##           Coding.region.change              Amino.acid.change
## 1 ENSCAFT00000000001:c.1802C>A ENSCAFT00000000001:p.Pro601His
## 2   ENSCAFT00000000003:c.12A>G   ENSCAFT00000000003:p.Ile4Met
##   Other.variants.within.codon Non.synonymous         annotation
## 1                          No            Yes          NOT FOUND
## 2                          No            Yes ENSCAFT00000000003
```

```r
 # produces an annotated CLC file with only the variants that were allocated
 # to a gene being retained.
AnnotCLCfile2 <- annot(a, genBED, type="BED", nomatch=FALSE, CLC=TRUE)
nrow(AnnotCLCfile2)
```

```
## [1] 9
```

```r
AnnotCLCfile2[c(1:2),]
```

```
##   Chromosome Region Type Reference Allele Reference.allele Length Linkage
## 2       chr1 363618  SNV         A      G               No      1      NA
## 3       chr1 604894  SNV         A      G               No      1      NA
##       Zygosity Count Coverage Frequency Probability Forward.read.count
## 2 Heterozygous    27       52    51.923           1                 27
## 3   Homozygous   212      212   100.000           1                134
##   Reverse.read.count Forward.reverse.balance Average.quality
## 2                  3                   0.100          36.067
## 3                 87                   0.394          36.204
##          Coding.region.change              Amino.acid.change
## 2  ENSCAFT00000000003:c.12A>G   ENSCAFT00000000003:p.Ile4Met
## 3 ENSCAFT00000000008:c.775T>C ENSCAFT00000000008:p.Phe259Leu
##   Other.variants.within.codon Non.synonymous         annotation
## 2                          No            Yes ENSCAFT00000000003
## 3                          No            Yes ENSCAFT00000000008
```
As demonstrated, this results in a difference of one row, due to the removal of the first variant in the second example. If nomatch=TRUE, variants that are not annotated, get the "NOT FOUND" flag in the annotation column. 

## Database filtering
Mendelian supports filtering against popular variant databases, like dbSNP. Three functions are available:
* a preparatory step: *prepvar* and *prepvarpar*
* the actual filtering: *varfilter*

*Prepvar* and *prepvarpar* require the input of a dbSNP file in a certain format. To make sure that the function works properly, the dbSNP file should be downloaded from the [UCSC website](http://genome.ucsc.edu/cgi-bin/hgTables?command=start), specifying the output format as "all fields from selected table" in the table browser. Take into consideration that downloading the complete dbSNP database can lead to large datafiles. The difference between *prepvar* and *prepvarpar* is that *prepvarpar* supports parallel computing. This is especially useful when the whole dbSNP is used. 

Important: before you read in the dbSNP file, you should remove the # in the header (before bin). If this is forgotten, an error will occur when reading in the file.

```r
dim(SNP)
```

```
## [1] 472  26
```
After preprocessing the vcf file or using the CLC file directly and loading the file containing the databasevariants, the databasevariantsfile is preprocessed:

```r
 y <- SNP
 reference <- "refNCBI"
 MAF <- 0.03
 a <- prepvar(y, MAF, reference) 
```

```
## 5 procent prepared 
## 10 procent prepared 
## 20 procent prepared 
## 30 procent prepared 
## 40 procent prepared 
## 50 procent prepared 
## 60 procent prepared 
## 70 procent prepared 
## 80 procent prepared 
## 90 procent prepared 
## Number of variants for filtering: 526
## Number of variants with > MAF: 355
```

```r
dim(a)
```

```
## [1] 355   4
```
Using a different MAF, the number of variants changes:

```r
 y <- SNP
 reference <- "refNCBI"
 MAF <- 0.7
 a <- prepvar(y, MAF, reference)
```

```
## 5 procent prepared 
## 10 procent prepared 
## 20 procent prepared 
## 30 procent prepared 
## 40 procent prepared 
## 50 procent prepared 
## 60 procent prepared 
## 70 procent prepared 
## 80 procent prepared 
## 90 procent prepared 
## Number of variants for filtering: 526
## Number of variants with > MAF: 73
```

```r
 dim(a)
```

```
## [1] 73  4
```
And finally, non-specifying MAF:

```r
 rm(MAF)
 y <- SNP
 reference <- "refNCBI"
 a <- prepvar(y,, reference) 
```

```
## 5 procent prepared 
## 10 procent prepared 
## 20 procent prepared 
## 30 procent prepared 
## 40 procent prepared 
## 50 procent prepared 
## 60 procent prepared 
## 70 procent prepared 
## 80 procent prepared 
## 90 procent prepared 
## Number of variants for filtering: 526
```

```r
 dim(a)
```

```
## [1] 526   4
```

After this preparatory step, the actual filtering is done:

```r
vcfvar <- VCFfile(test, "V10", TRUE, "PASS")
vcfvar
```

```
##      Chromosome  Region Allele     Zygosity
## 160        chr1 1563789      A   Homozygous
## 171        chr1 1581713      C Heterozygous
## 1711       chr1 1581713      G Heterozygous
## 176        chr1 1585597      G   Homozygous
## 177        chr1 1585642      T   Homozygous
## 178        chr1 1586752      C   Homozygous
```

```r
a <- prepvar(y, 0.03, reference) 
```

```
## 5 procent prepared 
## 10 procent prepared 
## 20 procent prepared 
## 30 procent prepared 
## 40 procent prepared 
## 50 procent prepared 
## 60 procent prepared 
## 70 procent prepared 
## 80 procent prepared 
## 90 procent prepared 
## Number of variants for filtering: 526
## Number of variants with > MAF: 355
```

```r
varfilter(vcfvar, a)
```

```
## Number of variants to be filtered: 6
## Number of variants from x after filtering retained: 1
```

```
##      Chromosome  Region Allele     Zygosity
## 1711       chr1 1581713      G Heterozygous
```

```r
a <- prepvar(y, 0.70, reference) 
```

```
## 5 procent prepared 
## 10 procent prepared 
## 20 procent prepared 
## 30 procent prepared 
## 40 procent prepared 
## 50 procent prepared 
## 60 procent prepared 
## 70 procent prepared 
## 80 procent prepared 
## 90 procent prepared 
## Number of variants for filtering: 526
## Number of variants with > MAF: 73
```

```r
varfilter(vcfvar, a)
```

```
## Number of variants to be filtered: 6
## Number of variants from x after filtering retained: 5
```

```
##      Chromosome  Region Allele     Zygosity
## 160        chr1 1563789      A   Homozygous
## 171        chr1 1581713      C Heterozygous
## 1711       chr1 1581713      G Heterozygous
## 176        chr1 1585597      G   Homozygous
## 177        chr1 1585642      T   Homozygous
```

```r
a <- prepvar(y,, reference) 
```

```
## 5 procent prepared 
## 10 procent prepared 
## 20 procent prepared 
## 30 procent prepared 
## 40 procent prepared 
## 50 procent prepared 
## 60 procent prepared 
## 70 procent prepared 
## 80 procent prepared 
## 90 procent prepared 
## Number of variants for filtering: 526
```

```r
varfilter(vcfvar, a)
```

```
## Number of variants to be filtered: 6
## No variants left from x after filtering
```

```
## [1] Chromosome Region     Allele     Zygosity  
## <0 rows> (or 0-length row.names)
```
The MAF used clearly affects the number of variants.

# Variant filtering
This section details on the main functions of the package: the variant filtering using different types of inheritance, penetrance and detectance. It can be broadly divided in two large sections: 
* nDom and nRec
* gDom and gRec

This division is based on the grouping unit that is required to be identical amongst cases. The first two functions look at the nucleotide level: all cases should have the same nucleotide (except under a reduced detectance). In the second group of functions, the cases are required to have an identical grouping unit (which is most often the same gene). The second group will always contain at least all the variants from the first group. A more detailed explanation can be found in the next sections. 

Before going into detail, a short explanation on penetrance and detectance, as it is used in these 4 functions:
* Detectance: 
    * mathematical: P(genotype|phenotype)
    * words: the probability of identifying a certain genotype, given the phenotype
    * calculation: the number of cases with the phenotype and the genotype, divided by the total number of cases.
    * It answers the question: if you have the disease, what is your genotype?
* Penetrance: 
    * mathematical: P(phenotype|genotype)
    * words: the probability of seeing a certain phenotype, given the genotype
    * calculation: the number of cases, divided by the sum of cases and the number of controls with the genotype but without the phenotype
    * It answers the question: if you have the mutation, what is the chance of you having the disease as well?
    
Important: the penetrance level is influenced by the detectance level chosen earlier. If the detectance is reduced (e.g genetic heterogeneity is present), the penetrance is calculated based on the number of individuals with a shared genotype.

## nDom and nRec

The first two functions discussed are nDom and nRec. As mentioned earlier, their major difference with the other two functions is the grouping unit. Here, the nucleotide is considered to be the unit of interest. 

For nDom, every non-reference variant (homozygous or heterozygous) present in a case can cause a disease. However, for nRec, only homozygous variants are able to cause a disease. This is the basic concept of Mendelian disorders. 

If the number of cases is > 1 and assuming 100% detectance, every case is assumed to have the same mutation, albeit homozygous for nRec and homo- or heterozygous for nDom. If the detectance drops however, not every person with a specific phenotype has necessarily the same mutation. It might be that the person is a phenocopy (= the organism's phenotype matches a phenotype which is determined by genetic factors, but in this case, it is enviromentaly induced) or that allelic heterogeneity (= the phenomenon in which different mutations at the same locus cause a similar phenotype) or locus heterogeneity (= the phenomenon in which mutations at different loci cause a similar phenotype) is present.

Adding controls increases the complexity. Under 100% penetrance, for nDom, every control is assumed to not having the disease causing mutation. This has the nice consequence that, for nDom, every variant in every control can be used for filtering. For nRec, only homozygous variants present in the controls can be used for filtering variants in the cases. If the penetrance drops beneath 100%, this means that some of the controls might have the mutation, but are not sick. There are numerous reasons for a reduced penetrance: age-dependent penetrance (the individual has to be old enough to get the disease, e.g. degenerative myelopathy in the dog), exercise-dependent penetrance (the individual has to be pushed far enough to show the phenotype, e.g. exercise-induced collapse in the dog), ... This makes filtering a lot more difficult and requires the mutation to be present in a sufficient number of controls before it can be used for filtering.

Finally, two family options are available for both nRec and nDom, adding additional assumptions. The two family options are:
* Ps-F: if both parents are available
* P-F: if only one of the two parents is available

In both cases, the parent(s) are assumed to be healthy  and the penetrance is assumed to be 100%. In addition, for nDom, the variant is assumed to be heterozygous in the progeny and not present in the parent(s). It is thus a reflection of the de novo mutation rate. For nRec, the variant is assumed to be homozygous in the progeny **and** heterozygous in the parent(s).

Some examples:
* nDom vs nRec with 2 cases:


```r
# restarting from the initial CLCfiles:
CLCfile1[,c(1:5,9)]
```

```
##    Chromosome             Region      Type Reference Allele     Zygosity
## 1        chr1             266523       SNV         G      T Heterozygous
## 2        chr1             363618       SNV         A      G Heterozygous
## 3        chr1             604894       SNV         A      G   Homozygous
## 4        chr1             699873       SNV         G      A   Homozygous
## 5        chr1      758487^758488 Insertion         -      G   Homozygous
## 22       chr1           13397637       SNV         T      A Heterozygous
## 23       chr1           13397669       SNV         C      T Heterozygous
## 48       chr1 16806440..16806442  Deletion       GAG      -   Homozygous
```

```r
CLCfile2[,c(1:5,9)]
```

```
##    Chromosome        Region      Type Reference Allele     Zygosity
## 1        chr1        266523       SNV         G      T Heterozygous
## 2        chr1        604894       SNV         A      G   Homozygous
## 3        chr1        699873       SNV         G      A   Homozygous
## 4        chr1 758461^758462 Insertion         -      G   Homozygous
## 5        chr1        758465       SNV         C      T   Homozygous
## 21       chr1      13397637       SNV         T      A Heterozygous
## 22       chr1      13397669       SNV         C      T Heterozygous
## 48       chr1      16009193       SNV         T      A   Homozygous
## 50       chr1      17094403       SNV         C      T   Homozygous
```



```r
output <- nDom(c("CLCfile1","CLCfile2"))

Number of cases: 2
Number of controls: 0 
Finished case 1 
Finished case 2 
Which detectance level is required? 
level (1): 50% (1/2)
level (2): 100% (2/2)
Choose level:2

output
  Chromosome   Region Allele Number of Samples
1       chr1 13397637      A                 2
2       chr1 13397669      T                 2
3       chr1   266523      T                 2
4       chr1   604894      G                 2
5       chr1   699873      A                 2

output <- nRec(c("CLCfile1","CLCfile2"))

Number of cases: 2
Number of controls: 0 
Finished case 1 
Finished case 2 
Which detectance level is required? 
level (1): 50% (1/2)
level (2): 100% (2/2)
Choose level:2

output
  Chromosome Region Allele Number of Samples
1       chr1 604894      G                 2
2       chr1 699873      A                 2
```
Where nDom retains 5 variants, nRec only retains 2. Under 100% detectance, nDom simply retains every shared variant. nRec retains only shared variants that are also homozygous. 
* nDom vs nRec for 1 case, 1 control:

```r
output <- nDom("CLCfile1", "CLCfile2"  )
Number of cases: 1
Number of controls: 1
Finished case 1 
Which detectance level is required? 
level (1): 100% (1/1)
Choose level:1
Finished control 1 
Which penetrance level is required? 
level (0): 100% (1/(1+0))
level (1): 50% (1/(1+1))
Choose level:0

output
  Chromosome             Region Allele Number of Samples
1       chr1 16806440..16806442      -                 1
2       chr1             363618      G                 1
3       chr1      758487^758488      G                 1


output <- nRec("CLCfile1", "CLCfile2")
Number of cases: 1
Number of controls: 1
Finished case 1 
Which detectance level is required? 
level (1): 100% (1/1)
Choose level:1
Finished control 1 
Which penetrance level is required? 
level (0): 100% (1/(1+0))
level (1): 50% (1/(1+1))
Choose level:0

output
  Chromosome             Region Allele Number of Samples
1       chr1 16806440..16806442      -                 1
2       chr1      758487^758488      G                 1
```
Where nDom retains 3 variants, nRec only retains 2. nDom removes every variant shared by case(s) and control(s). nRec only removes the variants that are homozygous in the controls. From the variants remaining in the case(s), it returns only the homozygous variants.

* The P-F option in nRec and nDom with 1 case, 1 control:

```r
output <- nDom("CLCfile1", "CLCfile2", "P-F")

Number of cases: 1
Number of controls: 1
P-F mode: 1 healthy parent, 1 affected progeny
Finished case 1 
Finished control 1 

output
  Chromosome             Region Allele Number of Samples
1       chr1 16806440..16806442      -                 1


output <- nRec("CLCfile1", "CLCfile2", "P-F")

Number of cases: 1
Number of controls: 1
P-F mode: 1 healthy parent, 1 affected progeny
Finished case 1 
Finished control 1 
No variants retained

output
NULL
```
nDom now returns only one variant, nRec no variants. nDom requires the variants in the case to be heterozygous after the removal of every single variant found in the control(s). One might argue that demanding the heterozygous state is not completely correct. It is true that two identical de novo mutations at the same place might occur. However, this is highly unlikely. If one wants to make sure that this possibility is taken into account, the nDom can be used but without the family option. For nRec, as no single variant is homozygous in the case (= progeny) AND heterozygous in the control (= parent), no variants are retained.

A wide variation of combinations are possible, especially if reduced penetrance and detectance are taken into account. As detailing every single option would lead us to far, this is omitted. 

## gDom and gRec
The next stop in variant filtering includes the following two functions: gDom and gRec. In this section, the grouping unit is not a nucleotide, but something bigger. This is rather vague as it is user dependent. For example, it might be an exon, a gene or even an entire chromosome. Most often however, it is a gene and this will be used further throughout this section.

For only one case, gDom is similar to nDom: every single variant is a possible disease causing candidate. For gRec however, for a gene to be retained it has to have at least one homozygous variant **or** be compound heterozygous (= the condition of having two heterogeneous recessive alleles at a particular locus that can cause genetic disease in a heterozygous state).

If the number of cases is > 1 and the detectance is 100%, gDom looks for every gene (= grouping unit) that has at least one mutation in every single case. Thus, in contrast with nDom, it does not require the same mutation to be present in every single case. For gRec, it first retains only those genes with at least one homozygous variant or compound heterozygosity. Next, from this reduced list of genes, it retains only those genes present in every single case. Hence, both functions allow for a certain degree of genetic heterogeneity to be present.

Controls can be added to filter the variants in the cases further. Again, we first assume 100% penetrance. For nDom, every variant present in every control can be used to filter the variants in the cases. After these variants have been filtered from the cases, the remaining variants are again grouped per gene and only those genes with at least one variant in every single patient are retained. For nRec, there is a two step filtering with controls. First, every homozygous variant present in every single control is used to filter the variants from the cases. Second, per control, all pairwise combinations from heterozygous variants within each gene are made. These are used to filter all pairwise combinations from variants present in the cases. Next, the same happens as mentioned in the previous paragraph for gRec.

Finally, reduced penetrance and detectance result in the most complex situation. For reduced penetrance, the situation is identical to nRec/nDom: the mutation has to be present sufficiently frequent in the controls before it can be used to filter variants from the cases. With reduced detectance, the level is not the mutation, but the gene (or any other functional unit chosen as defined by the user): the gene has to be sufficiently present in the cases to be retained, but it is not required anymore to be present in every single case.

Some practical examples:
* Before filtering, we first preprocess the CLC files. Again, this is only required for gDom and gRec and it should **not** be used prior nDom and nRec.

```r
CLCfile1proc <-CLCfile(CLCfile1, "Coding.region.change", TRUE, "; ")
CLCfile2proc <-CLCfile(CLCfile2, "Coding.region.change", TRUE, "; ")
CLCfile1proc[,c(1,2,4,5,9,22)]
```

```
##      Chromosome             Region Reference Allele     Zygosity
## 1          chr1             266523         G      T Heterozygous
## 2          chr1             363618         A      G Heterozygous
## 3          chr1             604894         A      G   Homozygous
## 4          chr1             699873         G      A   Homozygous
## 5          chr1      758487^758488         -      G   Homozygous
## 22         chr1           13397637         T      A Heterozygous
## 22.1       chr1           13397637         T      A Heterozygous
## 23         chr1           13397669         C      T Heterozygous
## 23.1       chr1           13397669         C      T Heterozygous
## 48         chr1 16806440..16806442       GAG      -   Homozygous
##                   group
## 1    ENSCAFT00000000001
## 2    ENSCAFT00000000003
## 3    ENSCAFT00000000008
## 4    ENSCAFT00000037436
## 5    ENSCAFT00000000011
## 22   ENSCAFT00000043665
## 22.1 ENSCAFT00000048679
## 23   ENSCAFT00000043665
## 23.1 ENSCAFT00000048679
## 48   ENSCAFT00000000149
```

```r
CLCfile2proc[,c(1,2,4,5,9,22)]
```

```
##      Chromosome        Region Reference Allele     Zygosity
## 1          chr1        266523         G      T Heterozygous
## 2          chr1        604894         A      G   Homozygous
## 3          chr1        699873         G      A   Homozygous
## 4          chr1 758461^758462         -      G   Homozygous
## 5          chr1        758465         C      T   Homozygous
## 21         chr1      13397637         T      A Heterozygous
## 21.1       chr1      13397637         T      A Heterozygous
## 22         chr1      13397669         C      T Heterozygous
## 22.1       chr1      13397669         C      T Heterozygous
## 48         chr1      16009193         T      A   Homozygous
## 50         chr1      17094403         C      T   Homozygous
##                   group
## 1    ENSCAFT00000000001
## 2    ENSCAFT00000000008
## 3    ENSCAFT00000037436
## 4    ENSCAFT00000000011
## 5    ENSCAFT00000000011
## 21   ENSCAFT00000043665
## 21.1 ENSCAFT00000048679
## 22   ENSCAFT00000043665
## 22.1 ENSCAFT00000048679
## 48   ENSCAFT00000000144
## 50   ENSCAFT00000000160
```
* gDom vs gRec with 2 cases:

```r
output <- gDom(c("CLCfile1proc","CLCfile2proc"))

Number of cases: 2
Finished case 1 
Finished case 2 
Which detectance level is required? 
level (1): 50% (1/2)
level (2): 100% (2/2)
Choose level:2
Number of variants retained: 10
Number of genes retained: 6

output
   Chromosome        Region Allele              Group Number of samples
1        chr1      13397637      A ENSCAFT00000043665                 2
2        chr1      13397637      A ENSCAFT00000048679                 2
3        chr1      13397669      T ENSCAFT00000043665                 2
4        chr1      13397669      T ENSCAFT00000048679                 2
5        chr1        266523      T ENSCAFT00000000001                 2
6        chr1        604894      G ENSCAFT00000000008                 2
7        chr1        699873      A ENSCAFT00000037436                 2
8        chr1 758461^758462      G ENSCAFT00000000011                 1
9        chr1        758465      T ENSCAFT00000000011                 1
10       chr1 758487^758488      G ENSCAFT00000000011                 1


output <- gRec(c("CLCfile1proc","CLCfile2proc"),,FALSE)

Number of cases: 2
Finished case 1 
Finished case 2 
Which detectance level is required? 
level (1): 50% (1/2)
level (2): 100% (2/2)
Choose level:2
Number of variants retained: 9
Number of genes retained: 5

output
  Chromosome        Region Allele              Group Number of samples
1       chr1      13397637      A ENSCAFT00000043665                 2
2       chr1      13397637      A ENSCAFT00000048679                 2
3       chr1      13397669      T ENSCAFT00000043665                 2
4       chr1      13397669      T ENSCAFT00000048679                 2
5       chr1        604894      G ENSCAFT00000000008                 2
6       chr1        699873      A ENSCAFT00000037436                 2
7       chr1 758461^758462      G ENSCAFT00000000011                 1
8       chr1        758465      T ENSCAFT00000000011                 1
9       chr1 758487^758488      G ENSCAFT00000000011                 1
```
These functions report both the number of functional units (genes here) and the number of variants retained. Here, gDom retains 10 variants in 6 genes, nRec only retains 9 variants in 5 genes.
* gDom vs gRec for 1 case, 1 control:

```r
output <- gDom("CLCfile1proc", "CLCfile2proc" )
Number of cases: 1
Number of controls: 1
Which detectance level is required? 
level (1): 100% (1/1)
Choose level:1
Finished control 1 
Which penetrance level is required? 
level (0): 100% (1/(1+0))
level (1): 50% (1/(1+1))
Choose level:0
Finished case 1 
Number of variants retained: 3
Number of genes retained: 3

output
  Chromosome             Region Allele              Group Number of samples
1       chr1 16806440..16806442      - ENSCAFT00000000149                 1
2       chr1             363618      G ENSCAFT00000000003                 1
3       chr1      758487^758488      G ENSCAFT00000000011                 1


output <- gRec("CLCfile1proc", "CLCfile2proc", FALSE )
Number of cases: 1
Number of controls: 1
Which detectance level is required? 
level (1): 100% (1/1)
Choose level:1
Finished control 1 
Which penetrance level is required? 
level (0): 100% (1/(1+0))
level (1): 50% (1/(1+1))
Choose level:0
Finished case 1 
Number of variants retained: 2
Number of genes retained: 2

output
  Chromosome             Region Allele              Group Number of samples
1       chr1 16806440..16806442      - ENSCAFT00000000149                 1
2       chr1      758487^758488      G ENSCAFT00000000011                 1
```
The same function with one control and one case returns 3 variants in 3 genes for gDom and 2 variants and 2 genes for gRec. 

## Special example: reduced detectance and reduced penetrance
Reduced penetrance and/or reduced detectance occur frequently, as discussed earlier. A practical example of how this comes in handy follows.

First, four different files are prepared from the same starting file.

```r
# cases: 2 cases with same variant, one with a different variant
case1 <-CLCfile1

a[,c(1:5,9)]
   Chromosome             Region      Type Reference Allele     Zygosity
1        chr1             266523       SNV         G      T Heterozygous
2        chr1             363618       SNV         A      G Heterozygous
3        chr1             604894       SNV         A      G   Homozygous
4        chr1             699873       SNV         G      A   Homozygous
5        chr1      758487^758488 Insertion         -      G   Homozygous
22       chr1           13397637       SNV         T      A Heterozygous
23       chr1           13397669       SNV         C      T Heterozygous
48       chr1 16806440..16806442  Deletion       GAG      -   Homozygous
case2 <-CLCfile1
case3 <-CLCfile1
case3[1,5] <- "C"
case3[,c(1:5,9)]
   Chromosome             Region      Type Reference Allele     Zygosity
1        chr1             266523       SNV         G      C Heterozygous
2        chr1             363618       SNV         A      G Heterozygous
3        chr1             604894       SNV         A      G   Homozygous
4        chr1             699873       SNV         G      A   Homozygous
5        chr1      758487^758488 Insertion         -      G   Homozygous
22       chr1           13397637       SNV         T      A Heterozygous
23       chr1           13397669       SNV         C      T Heterozygous
48       chr1 16806440..16806442  Deletion       GAG      -   Homozygous

# controls: 1 control
control1 <-CLCfile1
control1[1,5] <- "C"
```

Initially, only the cases are processed, including the case that has been changed (= case3). If we choose a detectance level of 100%, 7 of 8 variants are retained:


```r
nDom(c("case1","case2","case3"))
Number of cases: 3
Number of controls: 0 
Finished case 1 
Finished case 2 
Finished case 3 
Which detectance level is required? 
level (1): 33.3333333333333% (1/3)
level (2): 66.6666666666667% (2/3)
level (3): 100% (3/3)
Choose level:3
  Chromosome             Region Allele Number of Samples
1       chr1           13397637      A                 3
2       chr1           13397669      T                 3
3       chr1 16806440..16806442      -                 3
4       chr1             363618      G                 3
5       chr1             604894      G                 3
6       chr1             699873      A                 3
7       chr1      758487^758488      G                 3
```
If we allow the detectance to be reduced (level 2), all 8 variants are returned:

```r
nDom(c("case1","case2","case3"))
Number of cases: 3
Number of controls: 0 
Finished case 1 
Finished case 2 
Finished case 3 
Which detectance level is required? 
level (1): 33.3333333333333% (1/3)
level (2): 66.6666666666667% (2/3)
level (3): 100% (3/3)
Choose level:2
  Chromosome             Region Allele Number of Samples
1       chr1           13397637      A                 3
2       chr1           13397669      T                 3
3       chr1 16806440..16806442      -                 3
4       chr1             266523      T                 **2**
5       chr1             363618      G                 3
6       chr1             604894      G                 3
7       chr1             699873      A                 3
8       chr1      758487^758488      G                 3
```


# Combining filters
The final function was designed to further process the output of the four variant filter functions mentioned earlier. It assumes that all objects to be filtered, have the same phenotype (= input is assumed to be cases only, there is no room for controls) and that every object  has a column with an annotation. In addition, for the annotation to be useful, it should be the same for every single object: it does not make sense to group one object at the exon level and another object at the gene level. 

A practical example: 
* first: the preparation steps:

```r
output <- nRec("CLCfile1", "CLCfile2"  )
Number of cases: 1
Number of controls: 1
Finished case 1 
Which detectance level is required? 
level (1): 100% (1/1)
Choose level:1
Finished control 1 
Which penetrance level is required? 
level (0): 100% (1/(1+0))
level (1): 50% (1/(1+1))
Choose level:0

output
  Chromosome             Region Allele Number of Samples
1       chr1 16806440..16806442      -                 1
2       chr1      758487^758488      G                 1

output2 <- gDom("CLCfile1proc", "CLCfile2proc")
Number of cases: 1
Number of controls: 1
Which detectance level is required? 
level (1): 100% (1/1)
Choose level:1
Finished control 1 
Which penetrance level is required? 
level (0): 100% (1/(1+0))
level (1): 50% (1/(1+1))
Choose level:0
Finished case 1 
Number of variants retained: 3
Number of genes retained: 3

output2
  Chromosome             Region Allele              Group Number of samples
1       chr1 16806440..16806442      - ENSCAFT00000000149                 1
2       chr1             363618      G ENSCAFT00000000003                 1
3       chr1      758487^758488      G ENSCAFT00000000011                 1
```
We have created output for nRec and gDom now. Before we can continue processing, we have to annotate the output of the nRec function first:

```r
Annotoutput2 <- annot(output,genBED, type="BED", nomatch=FALSE, CLC=TRUE)

Annotoutput2
  Chromosome             Region Allele Number of Samples         annotation
1       chr1 16806440..16806442      -                 1 ENSCAFT00000000149
2       chr1      758487^758488      G                 1 ENSCAFT00000000011
```
* Now the objects are processed sufficiently for the *commonvar* function:

```r
x <- c("output2", "Annotoutput2")
group <- c("Group", "annotation")

a <-commonvar(x,group)
Which detectance level is required? 
level (1): 50% (1/2)
level (2): 100% (2/2)
Choose level:2
2 groups are retained 

a
  Chromosome             Region Allele              Group Freq
1       chr1 16806440..16806442      - ENSCAFT00000000149    2
2       chr1      758487^758488      G ENSCAFT00000000011    2
```

Some important remarks:
* zygosity is **not** important in this function: it assumes that that has been taken care of in the functions used before the *commonvar* function
* penetrance is not considered here as only cases are included in this function
* the detectance level is variable and works at the "functional unit" level (which is gene in this case).

# X-linked disorders
Altough not directly intented for this purpose, these functions can be used in specific cases when X-linked inheritance is expected:
* gDom: can be used every time an X-linked dominant inheritance is expected.
* nDom: idem, except for the the *family* option as with this option, the function expects the cases to be heterozygous and this is essentially only met when the case is female. If however the zygosity for a male case is heterozygous as well (by changing this manually or this might also be default depending on the program), the *family* option can be used as well for male cases. 
* nRec and gRec are more difficult as they make more assumptions towards zygosity. They can be safely used when all individuals (cases and if present, controls) are female. When males are included, two situations might occur:
    * The case is male: 
        * for nRec: mutations causing disease are assumed to be homozygous. As every single mutation can cause the disease in a male individual, changing the zygosity for every variant to homozygous (manually or this might also be default depending on the program), ensures this assumption to be met and a correct use of the function
        * for gRec: this is currently not supported. 
    * A control is male:
        * for nRec: this assumes again that only homozygous mutations can cause a disease. An identical operation as under nRec should solve this situtation: change the zygosity for the male control to homozygous for every single variant.  
        * for gRec: this is currently not supported.
        
In conclusion, the functions presented here can be used in almost every X-linked disorder, if the appropriate adjustments are made at the zygosity level for male individuals. 

# Overview of variant filters:

|filter | variant to be retained when shared by cases | variants in controls that will be used for filtering in cases |
|:-----:|:-------------------------------------------:|:-------------------------------------------------------------:|
|nDom | every variant | every variant |
|nRec | homozygous variants | homozygous variants |
|gDom | every gene with at least one variant | every variant |
|gRec | gene with at least one homozygous variant or at least once compound heterozygous | homozygous variants and/or pairwise combination of heterozgyous variants within each gene |

# Future plans:
* integration of Mendelian with other packages (e.g. VariantAnnotation)
* allow several separate analyses being run with one function call
* further increase the flexibility and options available in each function

# Flowcharts
## CLC file
![](http://www.heupdysplasie.ugent.be/Mendelian/figure3clc.png)

Remark: *annot* should be used between *nRec/nDom* and *commonvar*
## VCF file
![](http://www.heupdysplasie.ugent.be/Mendelian/figure4vcf.png)

Remark: *annot* should be used between *nRec/nDom* and *commonvar*
