Readme
================
yutian
2023-05-08

``` r
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("rtracklayer"))
BiocManager::install("rtracklayer")
if (!require("GenomicRanges"))
BiocManager::install("GenomicRanges")
library(rtracklayer)
library(dplyr)
library(GenomicRanges)
```

# preprocess

## Read JA.bed and Jb. bed

### For the first time

``` r
# set to the working directory
#if (!file.exists("JA.rds"))
#{JA <- read.table("../googledrivebedfiles/JA.bed",header = TRUE)
#colnames(JA) <- c("seqnames","start","end","beta","totreads","mreads","ureads")
#saveRDS(JA,"JA.rds")}
```

``` r
#if (!file.exists("JB.rds"))
#{JB <- read.table("../googledrivebedfiles/bed_files/JB.bed",header = TRUE)
#colnames(JB) <- c("seqnames","start","end","beta","totreads","mreads","ureads")
#saveRDS(JB,"JB.rds")}
```

### Every time after that

Since I already have the files, I will read them directly

``` r
# if (file.exists("JA.rds") & file.exists("JB.rds"))
# {JA <- readRDS("~/Desktop/DNAS/JASPAR/JA.rds")
# JB <- readRDS("~/Desktop/DNAS/JASPAR/JB.rds")}
```

## Filter the JA and JB files to only save rows with totreads \> 4.

``` r
# if (!file.exists("JAgreater4.rds"))
# {JAgreater4 <- subset(JA,JA$totreads>4)
# saveRDS(JAgreater4,"JAgreater4.rds")}
# if (!file.exists("JBgreater4.rds"))
# {JBgreater4 <- subset(JB,JB$totreads>4)
# saveRDS(JBgreater4,"JBgreater4.rds")}
```

Since those steps are already ran, and I have them saved, so I will just
read them directly.

``` r
JAgreater4 <- readRDS("~/Desktop/DNAS/JASPAR/JAgreater4.rds")
JBgreater4 <- readRDS("~/Desktop/DNAS/JASPAR/JBgreater4.rds")
```

## Match JA and JB ob the seqnames

``` r
sum(is.na(JAgreater4))
```

    ## [1] 0

``` r
sum(is.na(JBgreater4))
```

    ## [1] 0

``` r
JA_JB_inner_join <- inner_join(JAgreater4, JBgreater4, by = c("seqnames", "start"))
```

This will create a new data frame called JA_JB_inner_join that contains
only the rows that have matches between JA and JB based on the seqnames
and start columns. If there are any rows in JA or JB that do not have
matches in the other data frame, they will be excluded from the output.

``` r
sum(is.na(JA_JB_inner_join))
```

    ## [1] 0

Note that if there are any duplicate rows in either JA or JB based on
the seqnames and start columns, they will be included in the output
multiple times if they have matches in the other data frame.

``` r
# Remove duplicate rows from JA_JB_inner_join
JA_JB_inner_join <- distinct(JA_JB_inner_join, .keep_all = TRUE)
```

``` r
head(JA_JB_inner_join)
```

    ##   seqnames start end.x     beta.x totreads.x mreads.x ureads.x end.y    beta.y
    ## 1     chr1 13079 13080 0.80645161         31       25        6 13080 0.7500000
    ## 2     chr1 13216 13217 0.40000000         45       18       27 13217 0.4054054
    ## 3     chr1 13283 13284 0.83582090         67       56       11 13284 0.8947368
    ## 4     chr1 13302 13303 0.73972603         73       54       19 13303 0.8615385
    ## 5     chr1 13417 13418 0.04878049         41        2       39 13418 0.0000000
    ## 6     chr1 13503 13504 0.12500000         24        3       21 13504 0.2380952
    ##   totreads.y mreads.y ureads.y
    ## 1         24       18        6
    ## 2         37       15       22
    ## 3         57       51        6
    ## 4         65       56        9
    ## 5         32        0       32
    ## 6         21        5       16

## Read MA0003.2 from JASPAR as an example

``` r
# Load MA0003.2 bed file into GRanges object
ma0003.2_granges <- read.table("~/Desktop/DNAS/JASPAR/bed/MA0003.2.bed")
ma0003.2_granges <- GRanges(seqnames = ma0003.2_granges$V1, ranges = IRanges(start = ma0003.2_granges$V2, end = ma0003.2_granges$V3))
```

This will create a GRanges object called ma0003.2_granges that contains
the intervals from the MA0003.2 bed file.

``` r
head(ma0003.2_granges)
```

    ## GRanges object with 6 ranges and 0 metadata columns:
    ##       seqnames              ranges strand
    ##          <Rle>           <IRanges>  <Rle>
    ##   [1]     chr9   86577291-86577306      *
    ##   [2]     chr1 223168682-223168697      *
    ##   [3]     chr4     2050664-2050679      *
    ##   [4]    chr10 116394667-116394682      *
    ##   [5]    chr16   67982212-67982227      *
    ##   [6]    chr17   19967261-19967276      *
    ##   -------
    ##   seqinfo: 23 sequences from an unspecified genome; no seqlengths

``` r
# Resize intervals by adding 100 bp to each flank
ma0003.2_granges_resized <- resize(ma0003.2_granges, width = width(ma0003.2_granges) + 200, fix = "center")
head(ma0003.2_granges_resized)
```

    ## GRanges object with 6 ranges and 0 metadata columns:
    ##       seqnames              ranges strand
    ##          <Rle>           <IRanges>  <Rle>
    ##   [1]     chr9   86577191-86577406      *
    ##   [2]     chr1 223168582-223168797      *
    ##   [3]     chr4     2050564-2050779      *
    ##   [4]    chr10 116394567-116394782      *
    ##   [5]    chr16   67982112-67982327      *
    ##   [6]    chr17   19967161-19967376      *
    ##   -------
    ##   seqinfo: 23 sequences from an unspecified genome; no seqlengths

As you can see the ranges, there are 100 bps to each flank.

``` r
# Count number of intervals in MA0003.2
num_intervals_MA0003.2_resized <- length(ma0003.2_granges_resized)
num_intervals_MA0003.2_resized
```

    ## [1] 5098

As a result, there are 5098 intervals in the resized MA0003.2

``` r
# Load CpG bed file into GRanges object
JAB_granges <- GRanges(seqnames = JA_JB_inner_join$seqnames, ranges = IRanges(start = JA_JB_inner_join$start, end = JA_JB_inner_join$end.x))
head(JAB_granges)
```

    ## GRanges object with 6 ranges and 0 metadata columns:
    ##       seqnames      ranges strand
    ##          <Rle>   <IRanges>  <Rle>
    ##   [1]     chr1 13079-13080      *
    ##   [2]     chr1 13216-13217      *
    ##   [3]     chr1 13283-13284      *
    ##   [4]     chr1 13302-13303      *
    ##   [5]     chr1 13417-13418      *
    ##   [6]     chr1 13503-13504      *
    ##   -------
    ##   seqinfo: 388 sequences from an unspecified genome; no seqlengths

``` r
# Subset MA0003.2 intervals that overlap with CpG sites
ma0003.2_cpg_granges <- subsetByOverlaps(ma0003.2_granges_resized, JAB_granges)
# Count number of CpGs for each interval
num_cpgs<- countOverlaps(ma0003.2_granges_resized, JAB_granges)
head(num_cpgs)
```

    ## [1]  1  5  0  9 16  1

``` r
# Count number of intervals with 1, 2, 3, or more CpGs
num_intervals_1_cpg <- sum(num_cpgs == 1)
num_intervals_2_cpgs <- sum(num_cpgs == 2)
num_intervals_3_cpgs <- sum(num_cpgs == 3)
num_intervals_more_than_3_cpgs <- sum(num_cpgs > 3)
```

``` r
print(num_intervals_1_cpg )
```

    ## [1] 1172

``` r
print(num_intervals_2_cpgs)
```

    ## [1] 893

``` r
print(num_intervals_3_cpgs)
```

    ## [1] 515

``` r
print(num_intervals_more_than_3_cpgs)
```

    ## [1] 1210

There are 1173 intervals that have 1 CpG, 893 intervals that have 2
CpGs, 515 that have 3 CpGs and 1210 intervals that have more than 3
CpGs.

`ma0003.2_granges_resized` was created by resizing the original
ma0003.2_granges by adding 100bp to each end using the resize()
function.

The `subsetByOverlaps()` function was used to create
`ma0003.2_cpg_granges` by finding all the intervals in
`ma0003.2_granges_resized` that overlap with the intervals in
`JAB_granges`.

The `countOverlaps()` function was used to count the number of intervals
in `ma0003.2_granges_resized` that overlap with each interval in
`JAB_granges`. This creates a vector `num_cpgs` where each element
represents the number of overlapping intervals.

Finally, `num_intervals_1_cpg` is calculated as the sum of the elements
in `num_cpgs` that are equal to 1. This represents the number of
intervals in `ma0003.2_granges_resized` that overlap with exactly one
interval in `JAB_granges`.

``` r
# Find CpGs that overlap with MA0003.2 intervals
overlaps <- findOverlaps(ma0003.2_cpg_granges, ma0003.2_granges_resized)
# Get the interval names for each overlapping CpG
interval_names<- paste0(seqnames(ma0003.2_granges_resized)[subjectHits(overlaps)], "_", start(ma0003.2_granges_resized)[subjectHits(overlaps)], "..", end(ma0003.2_granges_resized)[subjectHits(overlaps)])
head(interval_names)
```

    ## [1] "chr9_86577191..86577406"    "chr1_223168582..223168797" 
    ## [3] "chr10_116394567..116394782" "chr16_67982112..67982327"  
    ## [5] "chr17_19967161..19967376"   "chr11_65342547..65342762"

``` r
names <- as.data.frame(table(interval_names))
head(names,20)
```

    ##               interval_names Freq
    ## 1    chr1_10003300..10003515    1
    ## 2  chr1_100128225..100128440    1
    ## 3  chr1_100147275..100147490    1
    ## 4  chr1_101468243..101468458    1
    ## 5  chr1_101723534..101723749    1
    ## 6    chr1_10488148..10488363    1
    ## 7  chr1_109104369..109104584    1
    ## 8  chr1_109817479..109817694    1
    ## 9  chr1_110334836..110335051    2
    ## 10 chr1_110334920..110335135    2
    ## 11 chr1_110878937..110879152    2
    ## 12 chr1_110878958..110879173    2
    ## 13 chr1_114302098..114302313    1
    ## 14   chr1_11544193..11544408    1
    ## 15 chr1_116520186..116520401    1
    ## 16 chr1_116652957..116653172    1
    ## 17     chr1_1167376..1167591    2
    ## 18     chr1_1167422..1167637    2
    ## 19 chr1_116832976..116833191    1
    ## 20 chr1_116876823..116877038    1

``` r
# Create a data frame to store the beta values for CpG sites within each interval
interval_beta <- data.frame(interval_names = interval_names, beta = numeric(length(interval_names)))
head(interval_beta)
```

    ##               interval_names beta
    ## 1    chr9_86577191..86577406    0
    ## 2  chr1_223168582..223168797    0
    ## 3 chr10_116394567..116394782    0
    ## 4   chr16_67982112..67982327    0
    ## 5   chr17_19967161..19967376    0
    ## 6   chr11_65342547..65342762    0

``` r
#write.csv(names, file="~/Desktop/DNA_Conservation/JASPARScript/intervanames_MA00032.csv",row.names = FALSE)
```

## Summary

This script is analyzing and processing genomic data.

After fileting the reads with greater tan 4 and stored them. The script
reads two `RDS files` (`JAgreater4.rds` and `JBgreater4.rds`) into two
data frames (`JAgreater4` and `JBgreater4`) and then performs an **inner
join** on the two data frames using the `inner_join()` function from the
`dplyr` package. The resulting data frame `(JA_JB_inner_join)` contains
only the rows that have matches between `JAgreater4` and `JBgreater4`
based on the seqnames and start columns.

The script then reads a BED file `(MA0003.2.bed`) into a `GRanges`
object `(ma0003.2_granges)` and resizes the intervals in the object by
adding 100 bp to each flank using the `resize()` function from the
`GenomicRanges` package. The resulting GRanges object
`(ma0003.2_granges_resized)` contains the intervals from the
MA0003.2.bed file with 100 bp added to each flank.

Next, the script creates another GRanges object `(JAB_granges)`which is
converted by `(JA_JB_inner_join)` by extracting the seqnames, start, and
end columns from `JA_JB_inner_join` and using them to create a new
GRanges object. The `subsetByOverlaps()` function from the GenomicRanges
package is then used to **subset ma0003.2_granges_resized to only
include intervals that overlap with the intervals in JAB_granges**. The
`countOverlaps()` function is used to count the number of intervals in
`ma0003.2_granges_resized` that overlap with each interval in
JAB_granges, and the resulting counts are used to calculate the number
of intervals in ma0003.2_granges_resized that have 1, 2, 3, or more
CpGs. Finally, the script prints the number of intervals with 1, 2, 3,
or more CpGs.

Overall, this script is analyzing the overlap between genomic intervals
and CpG sites, and using the GRanges and GenomicRanges packages to
manipulate and analyze genomic data.

## Table for a “paper”

### count the number of binding sites

``` r
num_sites_MA <- length(ma0003.2_granges_resized)
print(num_sites_MA)
```

    ## [1] 5098

``` r
num_sites_JAB <- length(JAB_granges)
print(num_sites_JAB)
```

    ## [1] 24882968

### number of TF binding sites with CpG in the interval

``` r
print(num_intervals_1_cpg )
```

    ## [1] 1172

``` r
print(num_intervals_2_cpgs)
```

    ## [1] 893

``` r
print(num_intervals_3_cpgs)
```

    ## [1] 515

``` r
print(num_intervals_more_than_3_cpgs)
```

    ## [1] 1210

The output shows the number of intervals with different numbers of CpGs
that overlap with the TF binding sites specified in the
`num_intervals_1_cpg`: the number of intervals that have only one CpG
overlapping with a TF binding site. `num_intervals_2_cpgs`: the number
of intervals that have two CpGs overlapping with a TF binding site.
`num_intervals_3_cpgs`: the number of intervals that have three CpGs
overlapping with a TF binding site. `num_intervals_more_than_3_cpgs`:
the number of intervals that have more than three CpGs overlapping with
a TF binding

### average methylation of the TF binding CpG (a histogram with X axis being percent methylation (0 to 100) and y axis being frequency)

``` r
# Subset MA0003.2 intervals that overlap with CpG sites
# Get the CpG locations in the interval
cpg_locs <- as.data.frame(ma0003.2_cpg_granges)
head(cpg_locs)
```

    ##   seqnames     start       end width strand
    ## 1     chr9  86577191  86577406   216      *
    ## 2     chr1 223168582 223168797   216      *
    ## 3    chr10 116394567 116394782   216      *
    ## 4    chr16  67982112  67982327   216      *
    ## 5    chr17  19967161  19967376   216      *
    ## 6    chr11  65342547  65342762   216      *
