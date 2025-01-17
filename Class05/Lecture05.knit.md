---
title: "Class05 Introduction to R Graphics"
author: "Liam Shweyk"
date: "January 24, 2019"
output: github_document
---


```r
# Lecture 5 R Graphics Intro

x <- rnorm(1000,0)
boxplot(x)
```

![](Lecture05_files/figure-markdown_github/unnamed-chunk-1-1.png)

```r
summary(x)
```

```
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## -3.06290 -0.65462  0.01326  0.01275  0.68623  3.58051
```

```r
hist(x)
```

![](Lecture05_files/figure-markdown_github/unnamed-chunk-1-2.png)

```r
boxplot(x, horizontal = TRUE)
```

![](Lecture05_files/figure-markdown_github/unnamed-chunk-1-3.png)

```r
# Hands on Section 2
weight <- read.table("bimm143_05_rstats/weight_chart.txt", header=TRUE)
plot(weight$Age, weight$Weight, typ="o", 
     pch=15, cex=1.5, lwd=2, ylim=c(2,10), 
     xlab="Age (months)", ylab="Weight (kg)", 
     main="Baby weight with age")
```

![](Lecture05_files/figure-markdown_github/unnamed-chunk-1-4.png)

```r
mouse <- read.table("bimm143_05_rstats/feature_counts.txt", sep="\t", header=TRUE)
barplot(mouse$Count, horiz = TRUE)
```

![](Lecture05_files/figure-markdown_github/unnamed-chunk-1-5.png)

```r
par(mar=c(3.1, 11.1, 4.1, 2))
barplot(mouse$Count, names.arg=mouse$Feature, 
        horiz=TRUE, ylab="", 
        main="Number of features in the mouse GRCm38 genome", 
        las=1, xlim=c(0,80000))
```

![](Lecture05_files/figure-markdown_github/unnamed-chunk-1-6.png)

```r
x <- c(rnorm(10000),rnorm(10000)+4)
hist(x, breaks=80)
```

![](Lecture05_files/figure-markdown_github/unnamed-chunk-1-7.png)

```r
mf <- read.delim("bimm143_05_rstats/male_female_counts.txt")
barplot(mf$Count, names.arg=mf$Sample, col=rainbow(nrow(mf)), 
        las=2, ylab="Counts")
```

![](Lecture05_files/figure-markdown_github/unnamed-chunk-1-8.png)

```r
barplot(mf$Count, names.arg=mf$Sample, col=c("blue2","red2"), 
        las=2, ylab="Counts")
```

![](Lecture05_files/figure-markdown_github/unnamed-chunk-1-9.png)

```r
genes <- read.delim("bimm143_05_rstats/up_down_expression.txt")
table(genes$State)
```

```
## 
##       down unchanging         up 
##         72       4997        127
```

```r
plot(genes$Condition1, genes$Condition2, col=genes$State,
     xlab = "Expression condition 1", ylab = "Expression condition 2")
```

![](Lecture05_files/figure-markdown_github/unnamed-chunk-1-10.png)

```r
palette(c("blue","gray","red"))
plot(genes$Condition1, genes$Condition2, col=genes$State, xlab="Expression condition 1", ylab="Expression condition 2")
```

![](Lecture05_files/figure-markdown_github/unnamed-chunk-1-11.png)

```r
# Let's plot expression vs gene regulation
meth <- read.delim("bimm143_05_rstats/expression_methylation.txt")
dcols <- densCols(meth$gene.meth, meth$expression)
plot(meth$gene.meth, meth$expression, col=dcols, pch=20)
```

![](Lecture05_files/figure-markdown_github/unnamed-chunk-1-12.png)

```r
# Find the indices of genes with above 0 expresion
inds <- meth$expression > 0

# Plot just these genes
plot(meth$gene.meth[inds], meth$expression[inds])
```

![](Lecture05_files/figure-markdown_github/unnamed-chunk-1-13.png)

```r
## Make a desnisty color vector for these genes and plot
dcols <- densCols(meth$gene.meth[inds], meth$expression[inds])

plot(meth$gene.meth[inds], meth$expression[inds], col = dcols, pch = 20)
```

![](Lecture05_files/figure-markdown_github/unnamed-chunk-1-14.png)

```r
dcols.custom <- densCols(meth$gene.meth[inds], meth$expression[inds],
                         colramp = colorRampPalette(c("blue2",
                                                      "green2",
                                                      "red2",
                                                      "yellow")) )

plot(meth$gene.meth[inds], meth$expression[inds], 
     col = dcols.custom, pch = 20)
```

![](Lecture05_files/figure-markdown_github/unnamed-chunk-1-15.png)

```r
plot(meth$promoter.meth, meth$gene.meth, ylab="gene methylation", xlab="promoter methylation")
```

![](Lecture05_files/figure-markdown_github/unnamed-chunk-1-16.png)

```r
# source the provided function so we can use it
source("bimm143_05_rstats/color_to_value_map.r")

mycols=map.colors(meth$expression, 
                  c(max(meth$expression), min(meth$expression)), 
                  colorRampPalette(c("blue","red"))(100))

plot(meth$promoter.meth, meth$gene.meth, 
     ylab="Gene Methylation", 
     xlab="Promoter Methylation", 
     col=mycols)
```

![](Lecture05_files/figure-markdown_github/unnamed-chunk-1-17.png)


---
title: "Lecture05.R"
author: "lshwe"
date: "Tue Feb 12 11:28:05 2019"
---
