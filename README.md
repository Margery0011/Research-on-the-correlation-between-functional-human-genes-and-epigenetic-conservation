README
================
Yutian Liu
2023-01-20

# Research on the correlation between functional human genes and epigenetic conservation

## Mentor

Name: Dr. Kimberly Siegmund

Email: <kims@usc.edu>

## Introduction

This project is based on the paper called [“Functional human genes
typically exhibit epigenetic
conservation”](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8439480/),
where they have found that DepMap critical genes exhibit cell-specific
preferential epigenetic conservation by comparing DNA methylation
measurements across intestinal crypts from the same person as well as
between duplicate cell lines.Therefore, they raised the hypothesis that
essential genes are more conserved because if there is a deviation, it
is not an important area of the genome since the body does not have to
protect it.

## Things have done

- Recauculated the relationship between Colon epithelial average gene
  expression (log2 CPM) and variability (variance) with gene
  conservation (PWD)
- Reproduced and optimize drawn figure for a better visualization of the
  relationship

``` r
library(readxl)
library(dplyr)
library(cowplot)
library(ggplot2)
library(ggpubr)
library(devtools)
library(Matrix)
library(ggridges)
library(data.table)
```

# Normal Colon data

Read Single cell RNA-seq data from normal colon (GSE125970 (38)).For
normal colon, all rectal and colon cells were analyzed, and all
epithelial cells were analyzed for the CRCs. Raw reads were converted to
CPM (reads divided by total reads) and then to log2CPM (x+1). Average
gene log2CPM expression and log2CPM variance were calculated for all
cells, including zero values.

Colon epithelial average gene expression (log2 CPM) and variability
(variance) can be calculated from single cell RNA-seq data.

Based on DepMap data \[1\] to identify “essential” genes.

``` r
nmx <- readxl::read_xlsx("Normal_colon.xlsx")
```

``` r
change_col_name <- function(df){
  colnames(df)[2] <-"Log2_CPM"
  colnames(df)[3] <-"Log2_CPM_variance"
  colnames(df)[4] <-"NN3_genePWD"
  return(df)
}
dfs <- list(Nona_NN)
Nona_NN <-lapply(dfs, change_col_name )
Nona_NN <- rbindlist(Nona_NN)
dfs <- list(Nona_essential_depmap)
Nona_essential_depmap <-lapply(dfs, change_col_name )
Nona_essential_depmap <- rbindlist(Nona_essential_depmap)
dfs <- list(Nona_NotDepMap)
Nona_NotDepMap <-lapply(dfs, change_col_name )
Nona_NotDepMap <- rbindlist(Nona_NotDepMap)
```

## DepMap Essential Genes

``` r
cor(Nona_essential_depmap$NN3_genePWD, Nona_essential_depmap$Log2_CPM, method = c("pearson", "kendall", "spearman"))
```

    ## [1] -0.2834464

``` r
cor.test(Nona_essential_depmap$NN3_genePWD, Nona_essential_depmap$Log2_CPM, method=c("pearson", "kendall", "spearman"))
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  Nona_essential_depmap$NN3_genePWD and Nona_essential_depmap$Log2_CPM
    ## t = -10.023, df = 1150, p-value < 2.2e-16
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.3357076 -0.2294457
    ## sample estimates:
    ##        cor 
    ## -0.2834464

If the correlation coefficient is greater than zero, it is a positive
relationship. Conversely, if the value is less than zero, it is a
negative relationship. A value of zero indicates that there is no
relationship between the two variables. From the result, negative and
significant correlations (p \< 0.05) were observed between gene PWDs and
gene expression/variability.

``` r
cor.test(Nona_essential_depmap$NN3_genePWD, Nona_essential_depmap$Log2_CPM, alternative = "less")
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  Nona_essential_depmap$NN3_genePWD and Nona_essential_depmap$Log2_CPM
    ## t = -10.023, df = 1150, p-value < 2.2e-16
    ## alternative hypothesis: true correlation is less than 0
    ## 95 percent confidence interval:
    ##  -1.0000000 -0.2382335
    ## sample estimates:
    ##        cor 
    ## -0.2834464

### Visualization

#### Figures Reproduced

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->
![](README_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

#### Density Plots & Change Y-X axis

![](README_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->
![](README_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

Many PWDs are concentrated in the range from 0.06 to 0.08

![](README_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

![](README_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

## Non-DepMap Essential Genes

``` r
cor(Nona_NotDepMap$NN3_genePWD, Nona_NotDepMap$Log2_CPM, method = c("pearson", "kendall", "spearman"))
```

    ## [1] -0.2655226

``` r
cor.test(Nona_NotDepMap$NN3_genePWD, Nona_NotDepMap$Log2_CPM, method=c("pearson", "kendall", "spearman"))
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  Nona_NotDepMap$NN3_genePWD and Nona_NotDepMap$Log2_CPM
    ## t = -31.747, df = 13288, p-value < 2.2e-16
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.2812547 -0.2496479
    ## sample estimates:
    ##        cor 
    ## -0.2655226

### Visualization

#### Figures Reproduced

![](README_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->
![](README_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

#### Density Plots & Change Y-X axis

![](README_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

![](README_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

Many PWDs are concentrated in the range from 0.05 to 0.15

![](README_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

![](README_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

## Comparison between Essential genes & Non-Essential genes

``` r
# Plot distributions
Nona_NN %>%
  ggplot(aes(x = Log2_CPM,
             y = `common essential DepMap`,
             fill = `common essential DepMap`)) +
  ggridges::geom_density_ridges(bandwidth = 4) + ggtitle("Essential, Non-essential gene expression comparison ")
```

![](README_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
# Compute the mean for and size of each group
group_means <-
  Nona_NN %>%
  group_by(`common essential DepMap`) %>%
  summarise(mean = mean(Log2_CPM),
            n = n())

group_means
```

    ## # A tibble: 2 × 3
    ##   `common essential DepMap`  mean     n
    ##   <chr>                     <dbl> <int>
    ## 1 DepMap_Essential           2.28  1152
    ## 2 Not_DepMap_Essential       1.16 13290

``` r
# Create our data visualisation
Nona_NN %>%
  ggplot(aes(x = `common essential DepMap`, y = Log2_CPM, fill = `common essential DepMap`)) +
  geom_boxplot() +

  # Add the mean for each group
  geom_point(data = group_means,
             aes(x = `common essential DepMap`, y = mean),
             shape = 3,
             size = 2)
```

![](README_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

``` r
# Compute the mean for and size of each group
group_means <-
  Nona_NN %>%
  group_by(`common essential DepMap`) %>%
  summarise(mean = mean(NN3_genePWD),
            n = n())
# Create our data visualisation
Nona_NN %>%
  ggplot(aes(x = `common essential DepMap`, y = NN3_genePWD, fill = `common essential DepMap`)) +
  geom_boxplot() +

  # Add the mean for each group
  geom_point(data = group_means,
             aes(x = `common essential DepMap`, y = mean),
             shape = 3,
             size = 2)
```

![](README_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

## Preliminary Conclusion

More variably and highly expressed genes are more conserved in normal
and neoplastic colon.

Colon epithelial average gene expression (log2 CPM) and variability
(variance) can be calculated from single cell RNA-seq data and
correlated (Pearson coefficient values in parenthesis) with gene
conservation (PWD). Negative and significant correlations (p \< 0.05)
were observed between gene PWDs and gene expression/variability.

PWDs in DepMap essential genes are concentrated in the range from 0.06
to 0.08 and PWDs in Non-DepMap essential gene are concentrated in the
range from 0.05 to 0.15.

## 1.20 Update Interactive Plots about TFs in paper

In the paper “Transcription factor expression as a predictor of colon
cancer prognosis: a machine learning practice”They identified five
transcription factors for the predictive model, and they are HOXC9,
ZNF556, HEYL, HOXC4 and HOXC6. Those five TFs(transcription factors) are
not in the common essential genes list. And in the paper “DNA
methylation events in transcription factors and gene expression changes
in colon cancer”, they selected 44 m-genes (methylation associated
genes) which showed a partial correlation adjusting by stromal content.
Those interactive plots are saved in this link. [Interactive
Plots](https://yutianl.shinyapps.io/DNAMethylationScatterPl/) You can
chose from the right legend and click the button twice to isolate them
for a better visualization. Choose the `Log2_CPM` as the Y variable and
`NN3_genePWD` as the X variable to view the association between gene
expression and pair distance.

## Things to do

- Download the transcription factor binding sites from [JASPAR
  database](https://jaspar.genereg.net/) to compare the regions we found
  interesting to examine if this hypothesis still stands in terms of
  transcriptome markers.
