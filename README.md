README
================
Yutian Liu
2023-02-20

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

-   Recauculated the relationship between Colon epithelial average gene
    expression (log2 CPM) and variability (variance) with gene
    conservation (PWD)
-   Reproduced and optimize drawn figure for a better visualization of
    the relationship

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
library(plotly)
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
expression and pair distance. 44 TFs in the paper are found in the colon
cancer data. Gene “RTEL1,” “PPP2R2C,” “WNT2,” “FAM150A,” “KRT23,” and
“FOXQ1” are not in the `Normal Colon` datasheet.

``` r
NN1<- readRDS("~/Desktop/DNA_Conservation/ScatterPl/census-app/data/NN1TFs220.rds")
plot_ly(data = NN1, x = ~NN3_genePWD, y = ~Log2_CPM,color = ~common.essential.DepMap,type = "scatter",
        mode   = 'markers') %>% 
        layout(title = "Normal Colon, see interactive plot at the website")
```

![](README_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

## 2.18 Update: Plot and Interactive Plots of ColonCancer data

These plots are drawn using the ColonCancer data.

### Interactive Plot

See it at the [website](https://yutianl.shinyapps.io/ColonCancer/)

On this website, the y-axis is `log2`, and you can choose the x variable
as you like. There are three categories, `DepMap` means those genes are
identified as more essential to proliferation in tissue culture in the
recent DepMap CRISPR-Cas9 single gene disruptions. `Not DepMap` means
those are not identified as essential DepMap genes. `PaperTFs` means
those listed as TFs in colon cancer publications. 50 TFs in the paper
are found in the `Colon Cancer` data.

``` r
coloncanShiny  <- readRDS("~/Desktop/DNA_Conservation/ColonCancer/ColonCancer/data/coloncanShiny3digits220.rds")
plot_ly(data = coloncanShiny, x = ~ColumnTinExcel, y = ~log2,color = ~DepMap,type = "scatter",
        mode   = 'markers') %>% 
        layout(title = "Colon Cancer, see interactive plot at the website")
```

![](README_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

These two figures show that TFs in the paper expressed more in colon
cancer than in the normal colon.

### With all data in the sheet

``` r
y <- coloncan$log2
x <- coloncan$ColumnTinExcel
ggplot(coloncan, aes(x = ColumnTinExcel, y = log2)) +
  geom_point() + 
  geom_smooth(method = "lm", colour = "green", fill = "green",formula = y~x) +
  theme_light() + stat_cor(method = "pearson", label.x = 0.10, label.y = 15)+ylim(-10,30) + ggtitle("ALL Colon Cancer, r =",round(cor(x,y),2))
```

![](README_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

As the plot shows, the Pearson correlation coefficient(r) is -0.55, and
the P value is less than 2.2e-16, indicating a significant result that
over all there is a negative correlation between log2 and Column T.

### With DepMap Data in Colon Cancer

``` r
print(cor.test(x = as.numeric(unlist(coloncanDep[20])), y = as.numeric(unlist(coloncanDep[,3]))))
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  as.numeric(unlist(coloncanDep[20])) and as.numeric(unlist(coloncanDep[, 3]))
    ## t = -1.1139, df = 1165, p-value = 0.2656
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.0898328  0.0248144
    ## sample estimates:
    ##         cor 
    ## -0.03261649

``` r
y <- coloncanDep$log2
x <- coloncanDep$ColumnTinExcel
ggplot(coloncanDep, aes(x = ColumnTinExcel, y = log2)) +
  geom_point() + 
  geom_smooth(method = "lm", colour = "green", fill = "green",formula = y~x) +
  theme_light() + stat_cor(method = "pearson", label.x = 0.05, label.y = 15)+ylim(-10,30) + ggtitle("Colon Cancer DepMap, r =",round(cor(x,y),2))
```

![](README_files/figure-gfm/unnamed-chunk-30-1.png)<!-- --> The plot
shows that the Pearson correlation coefficient(r) is -0.033. The P value
is less than 0.05, indicating a significant result that there is a
negative correlation between log2 and Column T in the DepMAP genes. This
correlation is less statistically significant than or as highly
correlated as those with all the data.

### With Not-DepMap Data in Colon Cancer

``` r
print(cor.test(x = as.numeric(unlist(coloncanNoDep[20])), y = as.numeric(unlist(coloncanNoDep[,3]))))
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  as.numeric(unlist(coloncanNoDep[20])) and as.numeric(unlist(coloncanNoDep[, 3]))
    ## t = -87.456, df = 18366, p-value < 2.2e-16
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.5523560 -0.5319349
    ## sample estimates:
    ##        cor 
    ## -0.5422255

``` r
y <- coloncanNoDep$log2
x <- coloncanNoDep$ColumnTinExcel
ggplot(coloncanNoDep, aes(x = ColumnTinExcel, y = log2)) +
  geom_point() + 
  geom_smooth(method = "lm", colour = "green", fill = "green",formula = y~x) +
  theme_light() + stat_cor(method = "pearson", label.x = 0.05, label.y = 20)+ylim(-10,25) + ggtitle("Colon Cancer NotDepMap, r =",round(cor(x,y),2))
```

![](README_files/figure-gfm/unnamed-chunk-33-1.png)<!-- --> This figure
is similar to the figure drawn with all of these colon cancer data,
there is a significant negative correlation between log2 and Column T,
and the r value is -0.54.

## 1.23 Update: Target gene candidate transcription factor prediction，HOXC4 as the example

JASPAR (<http://jaspar.genereg.net/>) is an open-access database
containing manually curated, non-redundant transcription factor (TF)
binding profiles for TFs across six taxonomic groups. They have released
the 9th version at 2022 which I am going to use.

### Obtain the base sequence of the potential promoter region of the target gene

Go to NCBI [gene database](https://www.ncbi.nlm.nih.gov/gene/), search
for HOXC4, then find the genomic location information of the gene in the
“Genomic context”directory. It shows that HOXC4 is located at Chr12:
54016888 - 54056030. ![Genomic
context](https://github.com/Margery0011/Capstone_Project/blob/main/Figures/fig1.png)

It is generally believed that the region 1000\~2000bp upstream of the
gene transcription start site is the promoter region of the gene. The
red arrow in the map below indicates that the gene is located in the
sense strand, and the direction of transcription is from left to right,
that is, the potential promoter region of the HOXC4 gene is
Chr12:54014888-54016888

Next, click the FASTA button in the Genomic regions, transcripts, and
products directory, enter the location information of the HOXC4 gene
promoter region on the right, and click Update View to get the potential
promoter sequence of the gene. Save the ![FASTA
result](https://github.com/Margery0011/Capstone_Project/blob/main/Figures/fig2.png)
for later.

### Predict potential transcription factors that bind to target gene promoter regions

Enter the [UCSC database](http://genome.ucsc.edu/) homepage, select
Track Hubs from the My Data drop-down menu and click to enter, enter
JASPAR in the Public Hubs search bar, click Search Public Hubs, find
Connect in the search results and click, the page refreshes to prompt
JASPAR Track loaded successfully.

Return to the UCSC homepage, select the genome version consistent with
the above from the Genomes drop-down menu, and the page is refreshed to
display the genome information browsing page. Click hide all to hide all
tracks, then select pack in the drop-down menu of the latest version of
JASPAR, and then click refresh in the upper right corner to set it to
only retain the information interface of JASPAR track.

Enter the position Chr12:54014888-554016888 of the potential promoter
region of the HOXC4 gene in the search box, and click GO to get the
potential transcription factors bound to the promoter region of the RET
gene. The direction of the arrow behind the transcription factor
indicates the transcription direction, and the transcription factor that
is consistent with the transcription direction of the target gene is
preferred; the color of the arrow behind the transcription factor
indicates the prediction score. The darker the color, the higher the
score, and the more reliable the corresponding prediction result. Click
JASPAR to enter the setting interface. Generally, the Minimum Score
above 200 can be regarded as statistically significant. Since there are
many predicted transcription factors, I set 500 as the Minimum Score
here, and the page refresh shows that the predicted transcription
factors are significantly reduced.

### Predict the binding site sequence of candidate transcription factors in the target gene promoter region

Take the four transcription factors ZNF148, ZNF384,MTF1whose arrows are
predicted to be darker in color as examples. ![UCSC Genome
Browser](https://github.com/Margery0011/Capstone_Project/blob/main/Figures/fig3.png)

Enter the JASPAR homepage, enter NFIC in the search box, set the same as
the above, select the latest version of the search result, tick the
front, and click Add to cart on the right. Add other transcription
factors to the shopping cart in the same way, and then click View cart.

On the right toolbar Scan, find the potential promoter region sequence
of the RET gene queried from the NCBI database above in FASTA format,
copy all of them and paste them into the search box. The threshold is
80% by default. When there are many prediction results, the threshold
can be increased. Set 85% here , click Scan, the page is refreshed, and
the result shows the predicted binding site sequence.

Prediction score, the higher the score, the more reliable the prediction
result; the start and end positions of the TFBS sequence; the specific
base sequence. Click Copy to save the prediction results.

## 0126 Update TFBSTools identifies transcription factor binding sites

### TFSTools interacts with JASPAR2022 to obtain PWM(position weight matrices)

Data needed :

PWM, position weight matrices, transcription factors binding sequence
(From Jaspar)

A DNA sequence, or fasta format sequence(From NCBI)

Here, I am using the 5 TFs(HOXC9, ZNF556, HEYL, HOXC4 and HOXC6)
mentioned before as an example to show the result.

``` r
opts <- list()
opts[["species"]] <- 'Homo sapiens'
opts["collection"] <-  'CORE'
PFMatrixList <- getMatrixSet(JASPAR2022, opts)
## To PWM
pwm <- toPWM(PFMatrixList)
pwm
```

    ## PWMatrixList of length 692
    ## names(692): MA0030.1 MA0031.1 MA0051.1 MA0059.1 ... MA1630.2 MA1633.2 MA0597.2

### Download sequence FASTA file from NCBI

### Gene Name Transfer

HOXC9, ZNF556, HEYL, HOXC4 and HOXC6 are symbol names, transfer them
into Entrezid for fetching them from NCBI.

``` r
E_genes <- c("HOXC9", "ZNF556", "HEYL", "HOXC4", "HOXC6")
gene_map <- select(org.Hs.eg.db, keys=E_genes, keytype="SYMBOL", columns=c("ENTREZID", "SYMBOL"))
```

    ## 'select()' returned 1:1 mapping between keys and columns

``` r
gene_map
```

    ##   SYMBOL ENTREZID
    ## 1  HOXC9     3225
    ## 2 ZNF556    80032
    ## 3   HEYL    26508
    ## 4  HOXC4     3221
    ## 5  HOXC6     3223

## 0130 update: Get the promoter sequence to find TFBS

``` r
transcriptCoordsByGene.GRangesList <-
     transcriptsBy (TxDb.Hsapiens.UCSC.hg19.knownGene, by = "gene") [gene_ids]
promoter.seqs <- getPromoterSeq (transcriptCoordsByGene.GRangesList,
                                   Hsapiens, upstream=2000, downstream=100)
```

``` r
gene_map
```

    ##   SYMBOL ENTREZID
    ## 1  HOXC9     3225
    ## 2 ZNF556    80032
    ## 3   HEYL    26508
    ## 4  HOXC4     3221
    ## 5  HOXC6     3223

``` r
HOXC4 <- promoter.seqs$`3221`
HOXC4sitesets <- searchSeq(pwm, HOXC4, seqname="HOXC4", min.score="80%", strand="+")
head(writeGFF3(HOXC4sitesets))
```

    ##                 seqname source feature start  end     score strand frame
    ## MA0030.1.3221      3221   TFBS    TFBS  1831 1844  9.033045      +     .
    ## MA0030.1.3221.1    3221   TFBS    TFBS   466  479  7.109460      +     .
    ## MA0030.1.3221.2    3221   TFBS    TFBS   477  490 11.821670      +     .
    ## MA0030.1.3221.3    3221   TFBS    TFBS   922  935  8.837639      +     .
    ## MA0031.1.3221.1    3221   TFBS    TFBS  1181 1188  8.475101      +     .
    ## MA0031.1.3221.2    3221   TFBS    TFBS  1497 1504  5.803724      +     .
    ##                                                                            attributes
    ## MA0030.1.3221   TF=FOXF2;class=Fork head/winged helix factors;sequence=ATTGGGTAAACATG
    ## MA0030.1.3221.1 TF=FOXF2;class=Fork head/winged helix factors;sequence=GTCTAATAAATAAA
    ## MA0030.1.3221.2 TF=FOXF2;class=Fork head/winged helix factors;sequence=AAAACATAAACAGT
    ## MA0030.1.3221.3 TF=FOXF2;class=Fork head/winged helix factors;sequence=TACATGTAAACACA
    ## MA0031.1.3221.1       TF=FOXD1;class=Fork head/winged helix factors;sequence=GTAAAAAG
    ## MA0031.1.3221.2       TF=FOXD1;class=Fork head/winged helix factors;sequence=CTAAAAAA

``` r
HOXC9 <- promoter.seqs$`3225`
HOXC9sitesets <- searchSeq(pwm, HOXC9, seqname="HOXC9", min.score="80%", strand="+")
head(writeGFF3(HOXC9sitesets))
```

    ##                 seqname source feature start  end    score strand frame
    ## MA0030.1.3225.1    3225   TFBS    TFBS  1184 1197 9.936481      +     .
    ## MA0030.1.3225.2    3225   TFBS    TFBS  1476 1489 6.720608      +     .
    ## MA0030.1.3225.3    3225   TFBS    TFBS  1944 1957 8.145674      +     .
    ## MA0031.1.3225.1    3225   TFBS    TFBS   938  945 5.991351      +     .
    ## MA0031.1.3225.2    3225   TFBS    TFBS   959  966 7.934533      +     .
    ## MA0031.1.3225.3    3225   TFBS    TFBS  1712 1719 5.803724      +     .
    ##                                                                            attributes
    ## MA0030.1.3225.1 TF=FOXF2;class=Fork head/winged helix factors;sequence=GAAAGATAAATAGA
    ## MA0030.1.3225.2 TF=FOXF2;class=Fork head/winged helix factors;sequence=CATACATACACACT
    ## MA0030.1.3225.3 TF=FOXF2;class=Fork head/winged helix factors;sequence=TTTGGGTAAATACG
    ## MA0031.1.3225.1       TF=FOXD1;class=Fork head/winged helix factors;sequence=ATAAAAAT
    ## MA0031.1.3225.2       TF=FOXD1;class=Fork head/winged helix factors;sequence=GTAAAAAC
    ## MA0031.1.3225.3       TF=FOXD1;class=Fork head/winged helix factors;sequence=TTAAAAAA

``` r
HOXC6 <- promoter.seqs$`3223`
HOXC6sitesets <- searchSeq(pwm, HOXC6, seqname="HOXC6", min.score="80%", strand="+")
head(writeGFF3(HOXC6sitesets))
```

    ##                 seqname source feature start  end     score strand frame
    ## MA0030.1.3223      3223   TFBS    TFBS  1831 1844  9.033045      +     .
    ## MA0030.1.3223.1    3223   TFBS    TFBS  1095 1108 13.688403      +     .
    ## MA0030.1.3223.2    3223   TFBS    TFBS  2063 2076  7.858104      +     .
    ## MA0031.1.3223.1    3223   TFBS    TFBS  1181 1188  8.475101      +     .
    ## MA0031.1.3223.2    3223   TFBS    TFBS  1497 1504  5.803724      +     .
    ## MA0031.1.3223.3    3223   TFBS    TFBS  1836 1843 13.755485      +     .
    ##                                                                            attributes
    ## MA0030.1.3223   TF=FOXF2;class=Fork head/winged helix factors;sequence=ATTGGGTAAACATG
    ## MA0030.1.3223.1 TF=FOXF2;class=Fork head/winged helix factors;sequence=CCAACATAAACAGG
    ## MA0030.1.3223.2 TF=FOXF2;class=Fork head/winged helix factors;sequence=CAGAAATAAATATT
    ## MA0031.1.3223.1       TF=FOXD1;class=Fork head/winged helix factors;sequence=GTAAAAAG
    ## MA0031.1.3223.2       TF=FOXD1;class=Fork head/winged helix factors;sequence=CTAAAAAA
    ## MA0031.1.3223.3       TF=FOXD1;class=Fork head/winged helix factors;sequence=GTAAACAT

``` r
HEYL <- promoter.seqs$`26508`
HEYLsitesets <- searchSeq(pwm, HEYL, seqname="HEYL", min.score="80%", strand="+")
head(writeGFF3(HEYLsitesets))
```

    ##                  seqname source feature start end    score strand frame
    ## MA0030.1.26508     26508   TFBS    TFBS   896 909 7.889733      +     .
    ## MA0031.1.26508.1   26508   TFBS    TFBS   245 252 5.991351      +     .
    ## MA0031.1.26508.2   26508   TFBS    TFBS   253 260 5.991351      +     .
    ## MA0031.1.26508.3   26508   TFBS    TFBS   593 600 5.803724      +     .
    ## MA0031.1.26508.4   26508   TFBS    TFBS   725 732 5.991351      +     .
    ## MA0031.1.26508.5   26508   TFBS    TFBS   885 892 5.803724      +     .
    ##                                                                             attributes
    ## MA0030.1.26508   TF=FOXF2;class=Fork head/winged helix factors;sequence=CCAAAACAAACAAA
    ## MA0031.1.26508.1       TF=FOXD1;class=Fork head/winged helix factors;sequence=ATAACCAT
    ## MA0031.1.26508.2       TF=FOXD1;class=Fork head/winged helix factors;sequence=ATAAATAT
    ## MA0031.1.26508.3       TF=FOXD1;class=Fork head/winged helix factors;sequence=TTAAAAAA
    ## MA0031.1.26508.4       TF=FOXD1;class=Fork head/winged helix factors;sequence=CTAAAAAT
    ## MA0031.1.26508.5       TF=FOXD1;class=Fork head/winged helix factors;sequence=CTAAAAAA

``` r
ZNF566 <- promoter.seqs$`80032`
ZNF566sitesets <- searchSeq(pwm, ZNF566, seqname="ZNF566", min.score="80%", strand="+")
head(writeGFF3(ZNF566sitesets))
```

    ##                  seqname source feature start end    score strand frame
    ## MA0030.1.80032     80032   TFBS    TFBS   828 841 6.323706      +     .
    ## MA0030.1.800321    80032   TFBS    TFBS   828 841 6.323706      +     .
    ## MA0031.1.80032.1   80032   TFBS    TFBS    38  45 8.475101      +     .
    ## MA0031.1.80032.2   80032   TFBS    TFBS    60  67 5.722194      +     .
    ## MA0031.1.80032.3   80032   TFBS    TFBS   107 114 5.991351      +     .
    ## MA0031.1.80032.4   80032   TFBS    TFBS   257 264 5.722194      +     .
    ##                                                                             attributes
    ## MA0030.1.80032   TF=FOXF2;class=Fork head/winged helix factors;sequence=CATAGGTACACAGA
    ## MA0030.1.800321  TF=FOXF2;class=Fork head/winged helix factors;sequence=CATAGGTACACAGA
    ## MA0031.1.80032.1       TF=FOXD1;class=Fork head/winged helix factors;sequence=GTAAGCAG
    ## MA0031.1.80032.2       TF=FOXD1;class=Fork head/winged helix factors;sequence=GTAAGAAA
    ## MA0031.1.80032.3       TF=FOXD1;class=Fork head/winged helix factors;sequence=TTAAAAAT
    ## MA0031.1.80032.4       TF=FOXD1;class=Fork head/winged helix factors;sequence=GTAAGAAA

## Things to do

-   Download the transcription factor binding sites from [JASPAR
    database](https://jaspar.genereg.net/) to compare the regions we
    found interesting to examine if this hypothesis still stands in
    terms of transcriptome markers.
