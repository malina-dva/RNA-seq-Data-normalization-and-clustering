
## Data quality assessment by sample clustering and visualization 

**For the first three steps of the tutorial (especially if you a newbie to R and RNA-seq data analyses), I recommend that you follow the YouTube Video in the link below, which was specifically created to demonstrate the installation of packages, setting working directory and downloading the input data set for this tutorial.**

[!["Getting started with RNA-seq"](http://img.youtube.com/vi/kR_iHVau8GI/0.jpg)]((https://www.youtube.com/watch?v=kR_iHVau8GI&t=1s))

You can find a thorough description of the input data set in step 3 in this tutorial.

Most of the code used here is available online:

https://bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#the-variance-stabilizing-transformation-and-the-rlog

I have obviously simplified everything a lot!

### 1. Installing and loading the DESeq2 package

First you need to load the package DESeq2, which is the **MAIN** package that will be used in this tutorial.
if you have never used DESeq2, you will have to install it beforehand.
The commands for installation you can find here:

https://bioconductor.org/packages/release/bioc/html/DESeq2.html

By the way this is how you can install not only DESeq2, but any package if you have R version >=3.6. Logically, you will have to replace DESeq2 with the name of the package you want to install. 

Alternatively you can install DESEq2 by typing: 

install.packages("DESeq2", dependencies = TRUE)

Regardless of the way you decide to install DESeq2, it will take some time for all the dependences and DESeq2 itself to be installed. So be patient. 

To load the package, so that you can make the functions in the package accessible, use the line below. 


```r
require("DESeq2")
```

**Lines starting with ## are an output of the code you ran earlier. No need to type or execute them**

```
## Loading required package: DESeq2
```

```
## Loading required package: S4Vectors
```

```
## Warning: package 'S4Vectors' was built under R version 3.6.1
```

```
## Loading required package: stats4
```

```
## Loading required package: BiocGenerics
```

You can also use an alternative command below to load the package, once installed:

library("DESeq2")

### 2. Setting working directory

Once we have DESEq2 installed and loaded, then we can set the working directory where you will have to
place all the files you will use as an input.
This is also the directory where all of the output files will be stored automatically once the directory is set.
Working directory can be any folder on your computer.

To set working directory in R studio, you use the function

setwd("")

You have to place the location of your working directory (folder) within the double quotes of this function.

Let's create a new **folder** (working directory) on your desktop and call it 

Tutorial_RNA_seq

To retrieve the full path to the location of the newly created folder, simply copy the whole folder and paste within the double quotes.
This is what I get when I do this:
"file:///C:/Users/Malina/Desktop/Tutorial_RNA_seq"
You will have to change the forward slashes to backward double slashes and you will have to remove everything before C: in order R to be able to use the path.
**NO NEED to change the type of slashes if you are using Linux and macOS machines.** 

Place the full path within the double quotes in the command below. 



```r
setwd("C:\\Users\\Malina\\Desktop\\Tutorial_RNA_seq")
```

### 3. Downloading and importing the RNA-seq data

**The RNA-seq data** which is going to be analyzed here is from this paper:

https://www.nature.com/articles/ncomms12418 

DNA hydroxymethylation controls cardiomyocyte gene expression in development and hypertrophy. 

The data is from cardiac myocytes (CMs) cells isolated from mice - embryonic day 14,5, neonatal, adult and TAC (transverse aortic constriction). The first three conditions represent developmental stages and the last condition is a hypertrophic (disease) state.

**Each of the four conditions has two biological replicates.**

In order to perform the analyses on this data set you will have to download the file with the raw counts (the input file) named Raw_counts_input.txt from github.
To do so go to:

https://github.com/malina-dva/RNA-seq-Data-normalization-and-clustering

Click on Raw_counts_input.txt file

Then click on Raw.

The file with the Raw_counts_input.txt will open in a separate page.

Right click and then click on "save page as" or "save as" to save the file, preferably on your Desktop. Make sure you save the file with its original name.

Once the file is saved on your Desktop move it to your working directory that we just created and set. (Tutorial_RNA_seq).


It is a good idea to open the input file in excel and check its content before you start to work with it in R. This is just for you to get an idea what type of format the data has. I don't recommend keeping the file open in excel while you are running the commands from the tutorial though.           
Repeating again! The file with the raw counts has to be placed in your working directory.
To read/load the input file in Rstudio use the command below.

```r
just.raw.counts = read.delim("Raw_counts_input.txt")
```

Let's check first several rows of the input file using the function 'head'.


```r
head(just.raw.counts)
```

```
##        Probe E14.5_R1 E14.5_R2 Neonatal_R1 Neonatal_R2 Adult_R1 Adult_R2 TAC_R1
## 1       Xkr4      229      363         545         417      133       96    280
## 2     Gm1992        0        3           0           0        0        0      1
## 3        Rp1        3        5           4          26       48       73      7
## 4      Sox17      206      195         285         226       51       44     53
## 5 AC129937.1        0        3           0           2        0        0      0
## 6     Mrpl15      597      599         468         480      318      316    389
##   TAC_R2
## 1    219
## 2      0
## 3     10
## 4     57
## 5      0
## 6    347
```


Now let's check the dimensions (total number of rows and columns) of your imported table. 


```r
dim(just.raw.counts)
```

```
## [1] 27195     9
```

There are 27195 genes (rows) and 9 columns. The first column (e.g. Probe) holds the names of the genes, the remaining 8 columns hold the gene counts (gene expression info) for the replicates of the four conditions we have. 
We want to actually specify that the column named 'Probe' has the information for the names of each row in the table. 

We can do that when we read the input file by specifying row.names =1, meaning that names of the rows are to be found in the first column.


```r
just.raw.counts = read.delim("Raw_counts_input.txt", row.names = 1)
```

```r
head(just.raw.counts)
```

```
##            E14.5_R1 E14.5_R2 Neonatal_R1 Neonatal_R2 Adult_R1 Adult_R2 TAC_R1
## Xkr4            229      363         545         417      133       96    280
## Gm1992            0        3           0           0        0        0      1
## Rp1               3        5           4          26       48       73      7
## Sox17           206      195         285         226       51       44     53
## AC129937.1        0        3           0           2        0        0      0
## Mrpl15          597      599         468         480      318      316    389
##            TAC_R2
## Xkr4          219
## Gm1992          0
## Rp1            10
## Sox17          57
## AC129937.1      0
## Mrpl15        347
```

After defining row.names = 1, the first column is no longer named Probe as you can see from output above.

```r
dim(just.raw.counts)
```

```
## [1] 27195     8
```
Also the number of columns is 8, not 9 - exactly matching the total number of replicates in the data set.

### 4. Generating the expression set for DESeq2

To do that we need three elements

1. names of the genes, which we already defined with row.names = 1

2. raw counts of the genes for the respective condition and replicate

it seems we already have the first two elements from the 'just.raw.counts' variable

3. metadata - phenotypic data for the expression set
The main purpose of the metadata is to define the relationship between replicates and conditions

The metadata is another text file called meta_data.txt.



Again download this file from github and place it in your working directory. 


To do so go to:

https://github.com/malina-dva/RNA-seq-Data-normalization-and-clustering

Click on meta_data.txt file

Then click on Raw

The file with the meta_data.txt will open in a separate page

Right click and then click on "save page as" or "save as" to save the file, preferably on your Desktop. Keep the original name of the file.

Move meta_data.txt to the Tutorial_RNA_seq folder.


Open it with excel to check its content.

Then let's load it in Rstudio.



```r
meta.data = read.delim(file="meta_data.txt", row.names = 1)
```
Again we specify the first column as row.names.


```r
head(meta.data)
```

```
##             condition
## E14.5_R1        E14.5
## E14.5_R2        E14.5
## Neonatal_R1  Neonatal
## Neonatal_R2  Neonatal
## Adult_R1        Adult
## Adult_R2        Adult
```
The column called 'condition' links the name of each replicate from 'just.raw.counts' to its respective condition.
For example E14.5_R1 and E14.5_R2, in the 'condition' column have the same name e.g E14.5, because they both belong to the E14.5 condition.


Now we can proceed to build DESeq object.

First we need to create DESeq matrix. 

```r
count.data.set <- DESeqDataSetFromMatrix(countData=just.raw.counts, 
                                         colData=meta.data, design= ~ condition) 
```

'countData' argument requires the table with the counts and the names of the genes

'colData' argument requires the metadata file

'design' requires the column of the metadata that holds the information for the relationships between replicates and conditions (e.g the 'condition' column from the metadata in our case)



Then we create DESeq object.

```r
count.data.set.object <- DESeq(count.data.set)
```

```
## estimating size factors
```

```
## estimating dispersions
```

```
## gene-wise dispersion estimates
```

```
## mean-dispersion relationship
```

```
## final dispersion estimates
```

```
## fitting model and testing
```

### 5. Data normalization and generation of dissimilarity tree

Having the DESeq object now we can check how replicates cluster based on the expression levels of the genes (counts)
using various clustering analyses.

We can use functions in R to calculate the degree of dissimilarity between replicates and conditions
based on the expression levels of the genes and then we can plot a tree of hierarchical clustering (dendrogram) between replicates.

Before we plot the dendrogram, it is a good idea
to normalize the data for 1. sequencing depth and 2. composition. 
Why do we need to do that?
1. Samples may have uneven amount of starting material. For example there may be samples that are difficult to obtain and the starting material for them would be much lower compared to the rest of the samples.
Furthermore some of the samples may be left on the sequencing machine longer than others. These factors may introduce unwanted variability in the sequencing counts between the compared samples which is not due to biological but technical reasons. We want to normalize for that sort of technical variability.  
2. If we are comparing samples from different tissues, we want to make sure that the genes expressed only in one or the other tissue are excluded when scaling factors are estimated (read below for more info). 

We will use 'vst' normalization (varianceStabilizingTransformation), which is part of the DESeq2 package.


```r
vsd <- vst(count.data.set.object)
```

Another way to normalize the data is using rlog function.

rld <- rlog(count.data.set.object)

You can read more about differences between vst and rlog normalization online at:

https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/varianceStabilizingTransformation

For now we will stick to vst normalization.
Once we have normalized the data we can export and save the transformed (normalized counts) to a text file.

To do that we first extract normalized counts in readable format from the DEseq object using the function "assay".


```r
norm.data = assay(vsd)
```

Now we can check the top few lines
of the newly created data frame (table) in the norm.data.

```r
head(norm.data)
```

```
##            E14.5_R1 E14.5_R2 Neonatal_R1 Neonatal_R2 Adult_R1 Adult_R2   TAC_R1
## Xkr4       7.677109 8.137802    8.915350    8.566363 7.849458 7.509247 8.505594
## Gm1992     4.330501 4.784299    4.330501    4.330501 4.330501 4.330501 4.675836
## Rp1        4.803272 4.914769    4.906902    5.752897 6.724628 7.200788 5.231634
## Sox17      7.552842 7.394836    8.080148    7.798572 6.784615 6.677513 6.606959
## AC129937.1 4.330501 4.784299    4.330501    4.739802 4.330501 4.330501 4.330501
## Mrpl15     8.886284 8.781576    8.714413    8.750078 8.962472 9.005624 8.936597
##              TAC_R2
## Xkr4       8.236562
## Gm1992     4.330501
## Rp1        5.418549
## Sox17      6.710847
## AC129937.1 4.330501
## Mrpl15     8.831688
```


You notice that the normalized counts for the genes are much smaller in comparison to the raw data.
This is due to fact that vst calculates scaling factor for each replicate considering seq depth and composition. All raw counts for the respective replicate are divided by the scaling factor for this replicate. 

When the scaling factors are determined, genes that have zero counts in the raw data in at least one, or more replicates are not considered. This is why vst normalization does not require raw data pre-filtering, by removing genes with low or zero counts.

VST has an upward shift for smaller values. In the current analyses zeroes are all replaced by the same value (e.g. 4.330501) in the norm.data.  


Then we save a table with the normalized data using the following command.


```r
write.table(norm.data, sep="\t",file="Norm_data_all_genes_NO_counts_cut_off.txt", row.names=TRUE,col.names=NA,quote=FALSE)
```



The newly created file will be in your working directory.

To perform hierarchical clustering analyses and to plot a dendrogram, we will evaluate dissimilarities (calculate Euclidean distance) between all eight replicates based on their normalized gene counts. 
We use the function below:


```r
sampleDists <- dist(t(norm.data),  method = "euclidean")
```

t in the command is 'transpose'. it reverses rows and columns in the table with the norm.data.

Let's have a look at what it does

```r
reversed_rows_columns = (t(norm.data))
```

Let's check first five rows and columns.

```r
reversed_rows_columns[1:5,1:5]
```

```
##                 Xkr4   Gm1992      Rp1    Sox17 AC129937.1
## E14.5_R1    7.677109 4.330501 4.803272 7.552842   4.330501
## E14.5_R2    8.137802 4.784299 4.914769 7.394836   4.784299
## Neonatal_R1 8.915350 4.330501 4.906902 8.080148   4.330501
## Neonatal_R2 8.566363 4.330501 5.752897 7.798572   4.739802
## Adult_R1    7.849458 4.330501 6.724628 6.784615   4.330501
```

To check the output with Euclidean distances we simply load the variable that we already generated earlier.


```r
sampleDists
```

```
##              E14.5_R1  E14.5_R2 Neonatal_R1 Neonatal_R2  Adult_R1  Adult_R2
## E14.5_R2     28.73316                                                      
## Neonatal_R1  84.10912  77.03313                                            
## Neonatal_R2  76.22477  70.24553    45.57274                                
## Adult_R1    149.68409 150.65134   128.50036   126.17128                    
## Adult_R2    154.88829 156.64020   136.58203   132.59935  36.11603          
## TAC_R1      135.37961 136.74983   115.48864   114.31138  52.83689  61.77608
## TAC_R2      137.99030 139.21268   115.15927   114.76331  47.06693  56.63482
##                TAC_R1
## E14.5_R2             
## Neonatal_R1          
## Neonatal_R2          
## Adult_R1             
## Adult_R2             
## TAC_R1               
## TAC_R2       34.88568
```

As you can see this is a matrix of the pairwise comparsions of the Euclidean distances for each replicate with any other.
The smaller the value, the smaller the difference between the replicates.
The replicates with the smallest difference are actually the E14.5 replicates.

Having the distance (dissimilarity) we can finally perform
hierarchical cluster analysis using hclust function.


```r
clusters=hclust(sampleDists)
```

Now we have everything to plot a dendrogram.

```r
plot(clusters)
```

![](https://github.com/malina-dva/RNA-seq-Data-normalization-and-clustering/blob/master/unnamed-chunk-22-1.png)<!-- -->

The replicates cluster together within their respective condition. Great!

### 6. Principal component analysis 

We can proceed to principal component analysis (PCA).
To explain properly the concept of PCA analysis I may need to creat a whole new tutorial. :)
But https://www.youtube.com/watch?v=HMOI_lkzW08 is a really good video regarding that!
plotPCA function (again part of the DESeq2 package) calculates and plots PCA for the first two principal components. We use the normalized counts from the vsd variable, similarly to the analyses above. Since we want replicates from the same condition to be colored identically on the PCA plot, we have to specify that the column "condition" from the metadata has to be taken as a information for replicates grouping. 



```r
plotPCA(vsd, intgroup=c("condition")) 
```

![](https://github.com/malina-dva/RNA-seq-Data-normalization-and-clustering/blob/master/unnamed-chunk-23-1.png)<!-- -->

Notice that in the legend, the conditions are ordered alphabetically, but we don't want that.
We want them to be organized by increasing developmental stage.

To do the latter we can use the package ggplot2 that has a function for reordering conditions in the legend.
First we load the package ggplot2.

If you don't have it installed you have to follow the same steps you used to install DESEq2, but replacing DESEq2 by ggplot2.


Once installed we load ggplot2


```r
require(ggplot2)
```

```
## Loading required package: ggplot2
```

```
## Warning: package 'ggplot2' was built under R version 3.6.2
```

Then we specify the desired order using scale_colour_hue function.


```r
plotPCA(vsd, intgroup=c("condition")) +
  scale_colour_hue(breaks = c("E14.5","E14.5","Neonatal","Neonatal","Adult","Adult","TAC","TAC"))
```

![](https://github.com/malina-dva/RNA-seq-Data-normalization-and-clustering/blob/master/unnamed-chunk-25-1.png)<!-- -->

Again the replicates cluster together in their respective condition and the trajectory of PC1 follows the trajectory of CMs development.

Importantly, the replicates for the TAC condition appear earlier than the replicates for Adult (healthy state) on the PC1 axis.
We know from other studies that the TAC condition activates neonatal gene program and thus it is logical that TAC appears closer to Neonatal than Adult on the PC1 axis.


Considering the current analyses we can conclude the following:

**1. There is no obvious technical variability in our samples (e.g. the data seems to be of a high quality)**

**2. Having performed these analyses we can proceed to identifying differentially expressed genes**























