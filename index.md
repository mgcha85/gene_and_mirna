# prediction target genes of miRNA using Lasso


Table of Contents
=================

  * [Data preparation](#Data-preparation)
  * [Pre-processing](#Pre-processing)
  * [Comparison](#Comparison)
  * [Methodology](#Methodology)
  * [Validation](#Validation)
  * [Point](#Point1)
  
  
  

Data-preparation
============
## Data-resources


|Resource | URL|
| ------------- | ------------- |
|Annotated genes | https://www.gencodegenes.org/human/release_32lift37.html|
|RNA-seq by tissues | https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-1733/samples/?s_page=1&s_pagesize=500|
|FANTOM by tissues | https://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz|

FANTOM CAGE peaks expression tables for the robust (expression>10TPM) set of peaks.  
It based expression table (RLE normalized) for human with annotation.  
RLE normalization is explained on [this paper](https://www.nature.com/articles/sdata2017112).  

## Converting Data Format

| Convert | Command |
| ------------- | ------------- |
| fastq → sam (hi-sat) | hisat2 -p 4 -x genome_tran -1 {fastq1} – 2 {fastq2} -S {sam} |
| sam → bam (samtools) | samtools sorted -@ 8 -o {bam} {sam} |
| bam → gtf (stringtie)	| stringtie -p 4 -G genes.gtf -o {gtf} -i {bam} |

## File entries by Tissues

|Resource | URL|
| ------------- | ------------- |
|FANTOM | https://drive.google.com/open?id=1y3w-rOvgbxtlCKWss1wAzcHD884Na__tijNWL7Yjb2I|
|RNA-seq | https://drive.google.com/open?id=1lKZN5hn2e6Zq3eKzn5r3H5NqdEfjgUmy|

## 22 common tissues between FANTOM and RNA-seq
The common tissues are listed in the below table. 


| RNA-seq | FANTOM |
| ------------- | ------------- |
| appendix | appendix |
| urinarybladder | bladder |
| brain	| brain|
| colon	| colon|
| esophagus | esophagus|
| fat | adipose|
| gallbladder | gall_bladder|
| heart | heart|
| kidney | kidney|
| liver | liver|
| lung | lung|
| lymphnode | lymph_node|
| ovary | ovary|
| pancreas | pancreas|
| placenta | placenta|
| prostate | prostate|
| salivarygland | salivary_gland|
| smallintestine | small_intestine|
| spleen | spleen|
| testis | testis|
| thyroid | thyroid|
| endometrium | uterus|


Pre-processing
============
## Averaged data by tissue
One tissue can have multiple replicates. Therefore, RNA-seq and FANTOM were averaged on same tissues. For this, we  only used chr1, chr2, …chr 22, chrX and chrY. The average is calculated by each chromosome and strand.

## Expression data by 22 tissues
RNA-seq data already contains FPKM while FANTOM does not have expression data. To get the expression data from FANTOM, we calculated expression level of (+/-)100 bp from annotated gene TSS. The expression level is summed by the score sum of cage tags within the region. When the summation, every tags, which are overlapped to the region;  [TSS-100, TSS+100] should be used.

## CAGE raw data vs. CAGE tag processed data
### Before  vs. After Filter

There are three filters.
1. maximum value of FANTOM or RNA-seq on 22 tissues is zero.  
2. median of FANTOM or RNA-seq on 22 tissues is zero.  
3. maximum value is greater than sum of others.  

### after filter
After the filtering, the total number of transcripts for 27,493 out of 84,088, which is the number of TSS in gencode.


## High consistent genes
RNA-seq and FANTOM data have expression level by a gene on the 22 tissues. Each gene has two 22 length vectors (RNA-seq and FANTOM), which an element show expression level on a tissue. By the two vectors, the correlation coefficient is calculated to observe how consistent between RNA-seq and FANTOM on a gene. Each gene has the coefficient and high coefficient shows two data source are consistent on the 22 tissues.  
For the correlation, **spearman method** was used. We round the coefficient at the decimal point with two digit.  
To extract only high correlated one, we set **threshold as 0.75**. Finally, **5,116 transcripts** are received as high consistent ones from **27,493** transcripts. the transcripts are grouped by gene name and if there are multiple transcripts, we picked one with the maximum score. After this, we got **2,312 genes**.  
The list of 2,312 transcripts is in the below link. only double type is used for every calculation. Transcript region are only used that gene and transcript type are protein coding gene from genecode transcripts.  
[link to see](https://drive.google.com/file/d/1Q9PvJdm1jVIW7zFH2ZE00rXkXMqLzZm4/view?usp=sharing)


## High consistent miRNA
For the consistent miRNAs, we compared eleven papers and retrieved the TSS locations by the miRNA name. The consistent miRNA TSSs are predicted by all papers and the TSS regions are close then we chose the meddle point as the TSS. The number of miRNA TSSs is 330. The name and TSS of 330 miRNAs are in the below link.  
[link to see](https://drive.google.com/open?id=1qZHBqubcYeJfVk7uvE8lFyzVXZeSknmIOf-qfS5zBmM)


Methodology
============
### FANTOM cell line data
we received cell line specific data. The below table shows the url of the data source and the table including file id, cell line and so on.

| label | URL|
| ------------- | ------------- |
| FANTOM url | http://fantom.gsc.riken.jp/5/datafiles/latest/basic/human.tissue.hCAGE/|
| cell line data list | https://drive.google.com/open?id=1frNwIbWzFwWdvKfThP8DT1mdbPfGRERy|

### Expression data by 240 cell lines
From the data preparation step, we got high consistent transcripts on gene and miRNA. We extracted the expression level from 240 cell lines. The score sum of the region where (+/-)100bp from annotated gene TSS is the expression level. We calculated this levels about the high consistent genes and miRNAs. 
some of 240 data does not have chromosomeY so, we did not consider the chromosome. Therefore, the total number of genes and miRNA is **2,313** and **330** as we mentioned. 

| label | URL|
| ------------- | ------------- |
| miRNA expression data url | https://drive.google.com/file/d/13pG6rZnK-Wn7Tn1cPrQe67-piAowJ5Jl/view?usp=sharing|
| gene expression data url | https://drive.google.com/file/d/12CPzPnH-d8bRTFvAtL-mUlsN6mQDnIJV/view?usp=sharing|


### Lasso linear regression
To predict target genes on a miRNA, we used lasso regression. the consistent genes and miRNAs have 240 vector, which are expression level on each cell line. Assume the number of genes and miRNA are m and n respectively. gene matrix **X**: (m x 240) and miRNA matrix **Y**: (n x 240).  
Lasso calculated the relationship between **X** and **Y**. The coefficient matrix by the Lasso result is **B** (n x m); **Y = BX**.  
Before, processing Lasso, **X** and **Y** were centered by subtracting mean value by each row.  

From the coefficient matrix (**B**), each row shows coefficient of one miRNA and multiple genes. Therefore, the non-zero coefficient values are assumed target genes, which are related to the corresponding miRNA. The coefficient of Lasso is [here](https://drive.google.com/file/d/1Y2sA1EE5KC4ZGyXqjyQnHA877p7NfC6O/view?usp=sharing)  

Figure 1 shows the mean of **B** as a sample.
![Image](/images/Figure_1.png)
Figure1. The mean of **B**    

When the coefficient is not zero between gene and miRNA From the result, we assume those are target genes by a miRNA.   
The target genes are [here](https://drive.google.com/file/d/1gQlfWrGh_GZEGHGTyeiVxajql9d7dSKL/view?usp=sharing).  
257 miRNAs have target genes out of 330. This is because the expression values on 240 cell lines from 73 miRNAs are all zeros.  
Therefore, 73 miRNAs failed to Lasso convergence.  

### Gene ontology analysis
From the lasso result, we got the target genes by a miRNA. To investigate the reliability of the result, we processed gene ontology analysis using the result. For the GO analysis, the high consistent genes are used as background genes. Target genes is the lasso result.  For the analysis, we use the below website.  
[link to see](http://cbl-gorilla.cs.technion.ac.il/)

From this analysis, we recevied another gene set based on the Lasso result.  
The result is [here](https://drive.google.com/file/d/1JXhQbQLv6scJV88k3499GNvHKP8zWBtk/view?usp=sharing).  

Figure 2 shows the input of GO by the website.
![Image](/images/Figure_2.png)  
Figure2. GO input set

![Image](/images/Figure_3.png)
Figure3. GO output  

### Gene Set Enrichment Analysis
For GSEA analysis, three inputs are required;  

1. Expression dataset  

2. phenotype labels  

3. gene sets  

For the expression dataset, we used the table which consist of gene rows and cell line columns. The gene names
of the dataset are required only Hugo symbol although GECODE contains
multi-symbols. Thus, the symbols were converted to Hugo symbols by
HGNC.  [https://www.genenames.org/tools/multi-symbol-checker/](https://www.genenames.org/tools/multi-symbol-checker/) 

For phenotype label, we selected
continuous file format as peak profile. The peak values are obtained by
averaging expression levels by cell lines and converted to integer value.  

Finally, molecular signature database (MSigDb) is used for gene set. We
had target genes by using Lasso. Through GSEA,
we searched statistically important gene set.  

In this experiment, we used some specific parameters. Since the expression dataset is already gene level data, we did not collapse dataset to gene symbols. GSEA ranks the genes and analyzes ranked list of genes. For the rank, various options exist. In this experiment, we used continuous phenotypes and Pearson, Cosine, Manhattan and Euclidean are available for the phenotypes. For this model, Manhattan is proper ranking method since it calculates similarity through distance with only one direction. The GSEA software is v3.0 from broadinstitute.  
The GSEA analysis result are in the below link.  

[link to see](https://drive.google.com/open?id=1SfNwJtYHWc1oRLbuy4lcu8Fd9cTsgpJj)  

Figure 4. shows input data of GSEA and Figure 5 shows the parameters of GSEA.  
![Image](/images/Figure_4.png)
Figure4. GSEA input profile  

![Image](/images/Figure_5.png)
Figure5. GSEA parameter set  

![Image](/images/Figure_6.png)
Figure6. GSEA heat map  


Validation
============
#### Lasso 
10 fold cross-validation  

For the cross-validation, we split the 240 vector to 216 (90%) as a train and 24 (10%) as a test. by the two set, we calculated the distance between actual Y and predicted Y. The below equation shows the distance.  

![Image](/images/formula1.gif)  
  
![Image](/images/formula2.gif)  

where Yactual is miRNA expression level, X is gene expression level and N is the size of Y matrix (n x 216 or n x 24). The below table shows 10-fold cross validation result.
(n: 330, N: 71,280 or 79,20)  


| # | test | train |
| ------------- | ------------- | ------------- |
| 0 | 0.0026 | 0.0230 |
| 1 | 0.0019 | 0.0232 |
| 2 | 0.0014 | 0.0236 |
| 3 | 0.0015 | 0.0237 |
| 4 | 0.0016 | 0.0232 |
| 5 | 0.0014 | 0.0238 |
| 6 | 0.0016 | 0.0234 |
| 7 | 0.0016 | 0.0237 |
| 8 | 0.0016 | 0.0236 |
| 9 | 0.0015 | 0.0235 |


### FANTOM cell lines vs. RNA-seq tissues
To validate 240 cell lines data, we checked out RNA-seq data as well. From EMBL-EBI, we have RNA-seq data on 27 tissues. The number of high consistent transcripts of gene and miRNA are m and n respectively. Then, X: (m x 27) and Y: (n x 27). The RNA-seq data contains only 264 miRNAs out of 330 miRNAs. By using this **X and Y**, we calculated coefficient matrix **B** by Lasso.  

To compare FANTOM and RNA-seq data, we calculated three distances as the below formula.

distance1: ![Image](/images/formula1-1.gif) and ![Image](/images/formula1-2.gif)  
distance2: ![Image](/images/formula2-1.gif) and ![Image](/images/formula2-2.gif)  

where B_{CELL} is calculated by using the 264 miRNAs’ expression data on FANTOM 240 cell lines, X_{RNA} is RNA-seq gene expression matrix with size (m x 27), Y_{RNA} is RNA-seq miRNA expression matrix with size (n x 27).  

 

### Lasso result vs. other researches
we use three researches; MIRANDA, RNA22 and TargerScan for the comparison. From the three researches, we extracted target genes by a miRNA. One version is the intersection of three researches and another version is the union of the researches. For example, miRNA1 has gene1, gene2 and gene3 by research1, gene2, gene4 and gene5 by research2 and gene2, gene6 and gene7 by research3. Then the intersection of the researches is gene2 and the union of the researches is gene1, gene2, gene3, gene4, gene5, gene6 and gene7. We applied GSEA and GO on this result.

We checked the correlation of expression level between RNA-seq and FANTOM tag using target genes of other researches. The union of target genes are used for this. The expression level are obtained as we have done in *"Expression data by 240 cell lines"*. 

We checked the Lasso regression if the target genes have many non-zeros values from the coefficient. Based on the result, only 4.06% target genes has non-zeros coefficient from intersection result and 3.55% target genes has non-zeros coefficient from union result

	
| | % non-zeros | mean of correlation |
| ------------- | ------------- | ------------- |
| Intersection | 4.06% | 0.02 |
| Union | 3.55% | 0.02 |


### Hypertest
we calculated hypergeometic test by comparing this result to other research or other result such as RNA-seq or GO. Since each miRNA has target genes, we can compare the target genes based on each miRNA. For this, we use survival function to calculate p-value. The function has four inputs, which q, m, n and k; q is # of the intersection of predicted genes by Lasso and predicted genes by another research. m is # of genes by another research, n is the # of high correlated genes – m and k is # of target genes by Lasso result. For the survival function, q-1 is used instead of q because the function is the inverse of the cumulative distribution function (1-cdf). Therefore the survival function excludes till q-1 instead of q.  

Not only we calculated p-value between this result and another result but we also consider the common genes by all three researches by a miRNA and common genes, which are overlapped by at least two methods and the union of the target gene are calculated. Also, GO has q-value by a miRNA so it is used for this comparison.  

The below link contains the all comparison for the hypergeometric test.  
[link to see](https://drive.google.com/file/d/1CA7FK5FSHs_Gd3LjONZLau5QnFNSILx_/view?usp=sharing)  


q significant is the multiplication of q-value (GO) and row number.  

The # significant column shows statistically importance of miRNA. The number shows how many research is important. The importance is considered if the p-value is smaller than 1 / # miRNA; threshold. Therefore, 2 means p-values of two researches are smaller than threshold and important.  
The result shows only 13 miRNAs are statistically important out of 327.  

Point1
=========

## Point 1-a
comparison of different TSS regions (500, 300, 100)  
Table to show for the three different TSS region sizes, the miRNA associated genes (from the lasso non-zero coefficients) are similar. For each miRNA, we compared between 500 and 300, 300 and 100, 500 and 100, how many percent of their associated genes are the same? For each percentage, you will use the number of shared genes to divide the smaller number of genes for a region size, In 300 case. we made a table for this so that in each row, you have the miRNA name, the three percentages for the three comparisons.

To compare the different size of mirRNA TSS region, we got target genes by Lasso regression on (+/-)100, 300, 500 respectively.  
#100, 300, 500 means the number of size 100, 300, 500 target genes. ∩/∪ means intersection genes of the pair (eg: 100 vs. 300) divided by union genes of the pair. 
∩/# smaller means intersection genes of the pair and divided by the smaller number of either of the pair.  
The result is [here](https://drive.google.com/file/d/1kWHfbK46wwoB-L5tmn4umf_kIq_zWy0f/view?usp=sharing).  

### statistics of a column with the larger number of the two numbers in each row of the three tables
(q20, q40, q60 is quantile 20%, 40%, 60%)
| pair | mean | median | min | max | q20 | q40 | q60|
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | 
| 100 vs 300 | 0.9639 | 1 | 0.6313 | 1 | 0.9427 | 1 | 1|
|300 vs 500| 0.9815 | 1 | 0.6537 | 1 | 0.9952 | 1 | 1|
|100 vs 500| 0.9513 | 1 | 0.6313 | 1 | 0.9107 | 1 | 1|


## Point 1-b

### summary of lasso regression result
The target genes from Lasoo regression are obtained by non-zero values.
The result contains miRNA name, target genes and the number of genes [here](https://drive.google.com/file/d/17FG5ArfkgSZoPBfT0ptUKUIfhn__W1X4/view?usp=sharing).

### statstics of the number of target genes
| mean | median | min | max | q70 | q80 | q90|
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | 
| 409.20 | 381 | 121 | 1282 | 439 | 484 | 606|


## Point 1-c

### Identified genes correlate with the miRNAs better than target genes
We compared Lasso result to two other researches (TargetScan, MirTarBase).
Each research has their own result. Each miRNA has target genes.

First, we converted their miRNA to pre-miRNA because we used consistent pre-miRNA TSS.
Second, gene name is also converted to transcript ID if they only provide gene name. In this process, we use all transcripts of corresponding gene.

After conversion, we calculated expression level by 240 cell-lines as we did in "Expression data by 240 cell lines".
Now, we have each expression values of miRNAs and genes.
A miRNA and gene have 240 vector, which is each cell line expression level.

From the expression level, we calculated correlation coefficient between miRNA and a target genes from each research.  
For example, a miRNA has 80 target genes,the miRNA has 80 length vector, which is correlation coefficient of the miRNA and each target gene.  
For correlation coefficient, we used spearman correlation.  

To compare the correlation coefficients among three research, we used Mann–Whitney U test.  
This test shows difference between the distributions of the data samples.  

The result is [here](https://drive.google.com/file/d/1RkeB0SrE29iVR-MSnlQZjw-QJDmXzFP_/view?usp=sharing).  
The table contains the following columns [stat (ts), p-value (ts), # genes (ts), stat (mi), p-value (mi), # genes (mi), # genes (Lasso)]

- stat (ts): U-test of TargetScan  
- p-value (ts): p-value of TargetScan  
- #genes (ts): the number of target genes from TargetScan  
- stat (mi): U-test of miRTarBase  
- p-value (mi): p-value of miRTarBase  
- #genes (mi): the number of target genes from miRTarBase  
- #genes (Lasso): the number of target genes from Lasso result.  

  
  
The below table shows how many miRNAs has smaller than p-value by two researches.
| research | # p<0.01 | # p<0.001 | # p<0.01 (genes>=10) | # p<0.001 (genes>=10) |
| ---- | ---- | ---- | ---- | ---- |
| TargetScan | 30 | 26 | 30 | 26 |
| miTarBase | 167 | 144 | 159 | 139 |

- #p<0.01: the number of miRNAs with p-value<0.01  
- #p<0.001: the number of miRNAs with p-value<0.001
- #p<0.01 (genes>=10): the number of miRNAs with p-value<0.01 and the number of target genes is at least 10.
- #p<0.001 (genes>=10): the number of miRNAs with p-value<0.001 and the number of target genes is at least 10.


## Point 1-d
### Identified genes are similar using tissue or cell line samples
We compared target genes using 93 tissues and 240 cell lines data, which is in section **"File entries by Tissues"**.
We checked how many target genes are shared by a miRNA. 
The result is in [here](https://drive.google.com/file/d/1rH7RL2rJxF7DODG89PSZ5QZAWcoTKR1E/view?usp=sharing).  

#### The column description
- #tissue: the number of target genes by 93 tissues  
- #cell_lines: the number of target genes by 240 cell lines  
- #∩: the number of intersection of two target genes
- #U: the number of union of two target genes  
- #∩/#U: the ratio of intersection to union  
- #∩/#small: the ratio of intersection to smaller number of #tissue and #cell_lines.
- #∩/#large: the ratio of intersection to bigger number of #tissue and #cell_lines.
- larger: larger one between #∩/#small and #∩/#large (in this case, #∩/#small is always larger)


The below table shows the statistics of larger  
| mean | median | min | max | q70 | q80 | q90 | 
| ---- | ----- | --- | --- | --- | --- | --- | 
| 0.5917 | 0.5864 | 0.4327 | 0.9196 | 0.6212 | 0.6478 | 0.6731 |



Point2
=========

## Point2-a. some have GO significance
There are three ontology mode "biological process, molecular function, cellular components".
We processed all modes and compared.  

[cellular component](https://drive.google.com/file/d/1G_aoYpLlh7rzcWHPbwBazO4h2xfKgAve/view?usp=sharing)  
[molecular function](https://drive.google.com/file/d/1iXdJEZMNka1IPPpf2XgyOOP7_hw_Q_g4/view?usp=sharing)  
[biological process](https://drive.google.com/file/d/17pYj7d4mlzy-r4gMeWyM25ADZvBR3Sqk/view?usp=sharing)  

And then we chceck their significance. for this, we checked corrected p-value <0.01 by each miRNA.  
If one of three at least has >0.01, the miRNA is not significant.  
The result if [here](https://drive.google.com/file/d/1JXhQbQLv6scJV88k3499GNvHKP8zWBtk/view?usp=sharing)  

The below table shows how many significant miRNA is by each mode out of 258 miRNA.  
| | cellular component | molecular function | biological process | 
| --- | --- | --- | --- |
| **#** | 254 | 217 | 244 |



## Point2.b. some have the target significance
I made one more table how many of the lasso genes are from other researches (targetscan and miTarBase) and how many percent of the lasso genes are from them.

The result is [here](https://drive.google.com/file/d/1DBqGdCx-RZR_ytQ4k5vrCgdaAmVtWBwd/view?usp=sharing)  

The below table shows the statistics of the above table.  
- mean of #: mean of the number of shared genes  
- mean of %: mean of the percent of shared genes  
and so on.

| research | mean of # | mean of % | median of # | median of % |  min of # | min of % | max of # | max of % | q70 of # | q70 of % | q80 of # | q80 of % | q90 of # | q90 of % | 
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **miTarBase** | 5.829 | 0.0143 | 4 | 0.0088| 1 | 0.0014 | 153 | 0.1839 | 5 | 0.0147 | 7 | 0.0176 | 10 | 0.025 |
| **TargetScan** | 13.55 | 0.034 | 11.5 | 0.026 | 1 | 0.0025 | 47 | 0.1 | 13 | 0.036 | 15 | 0.052 | 29 | 0.059 |


## Point2.c. compare with the GSEA gene sets.
- n: target genes by Lasso
- N: consistent genes
- m: Lasso genes ∩ GSEA gene set
- M: consistent genes ∩ GSEA gene set
- K: the number of miRNAs (257)

We calculated the p-value as 1-phyper(m-1, M, N-M, n).  
If this p-value is smaller than 0.01/K, we think the lasso genes of this miRNA is significantly overlapping with a GSEA gene set.  
we repeated this for every GSEA gene set (h, c1, ..., c8) and record those significant GSEA gene set names and their p-values. 
I got a table with the first column as the miRNAs, the second column is signficant and third column is p-values and phyper input parmameter in the fourth column.  
#genes (mi) and #genes (ts) mean >=10 common genes between mirtarbase/targetScan and the lasso genes of corresponding miRNA.  
The corresponding GSEA gene set name is sheet name. 

The table is [here](https://drive.google.com/file/d/1jVnaqS3lHda-QBkhXRj26Wa6oggOD93T/view?usp=sharing)  

The below table shows how many miRNAs are significant out of K by a gene set. 
| - | h | c1 | c2 | c3 | c4 | c5 | c6 | c7 | c8 | 
| -- | -- | -- | -- | -- | -- | -- | -- | -- | -- |
| **#** | 257 | 0 | 94 | 14 | 257 | 253 | 257 | 6 | 257 |

## Point2.d. overall of the miRNAs with at least one type of supporting evidence

[This table](https://drive.google.com/file/d/1Y99MYPKgaiHuvtA-32fcizAFdrGPA9nz/view?usp=sharing) is the summary for how many miRNAs, we have either significant GO or significant all types of GSEA. 

The table also include a miRNA has how many target genes from targetScan (#genes(ts)) or mirTarbase (#genes(mi)).  


## Point3.a. 10-fold cross-validation on the 240 cell lines

We trained the model with the 240x0.9=216 cell lines and then test the trained model on the remaining 24 cell lines.  
For each of the 10 tests, we used the trained matrix B from the 216 cell lines to the data in the remaining 24 cell lines to calculate the difference between the predicted miRNA expression and the actual miRNA expression in the 24 cell lines.   
  
![Image](/images/formula6.gif)    
![Image](/images/formula7.gif)    

In [this table](https://drive.google.com/file/d/1_7teMWi_IHfIjL2sKeqLPGiFuRJJatFk/view?usp=sharing), each cell shows  Y_{trn}, Y_{pred}, Y_{diff}.  
sheet name means cross validation number.  
- X_{tst): (2312 x 24)
- B_{trn}: (257 x 2312)
- Y_{pred}: (257 x 24)
- Y_{tst}: (257 x 24)

**trn: train, tst: test**  

The three last column in the table,  
- **avg(diff)**: avgerage(|Y_{tst} - Y_{pred}|)
- **avg(distance)**: the average of a pair from 24 cell lines. In this case, Y_{tst} has 24 cell lines.  

pair1: |cell line1 - cell line2|  
pair2: |cell line1 - cell line3|  
pair3: |cell line1 - cell line4|  
.  
.  
.  
pair276: |cell line23 - cell line24|  
The avg(distance) is the average of the 276 values.  

we count how many avg(diff) < avg(distance) of miRNAs are there.

The below table shows the result across 10 test sets.
|  | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
| - | - | - | - | - | - | - | - | - | - | - |
| # | 75 | 89 | 80 | 87 | 73 | 189 | 94 | 93 | 81 | 84 |

<br>
Also,  
The avg(distance) is |Y_{trn} - Y_{pred}| in the 216 training cell lines instead of the 24 testing cell lines.  
 
![Image](/images/formula6-1.gif)    
![Image](/images/formula7-1.gif)    

In this [table](https://drive.google.com/file/d/1MAW7xrooPExHp9FYR8zLtMwN18_IkOry/view?usp=sharing), each cell shows Y_{trn}, Y_{pred}, Y_{diff} as the above table.

The below table shows the result across 10 train sets.
|  | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
| - | - | - | - | - | - | - | - | - | - | - |
| # | 100 | 112 | 121 | 113 | 90 | 102 | 100 | 114 | 101 | 98 |


Another table is [here](https://drive.google.com/file/d/1erA1fyT0CLxiaLDLYxtnGTiYXBFcS8v3/view?usp=sharing).  
The ratio of average distance; Y_{ratio} = |Y_{pred} - Y_{tst}| / Y_{tst} is shown in the table. the format is same as the above tables.  
|  | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
| - | - | - | - | - | - | - | - | - | - | - |
| # | 83 | 105 | 93 | 81 | 107 | 86 | 87 | 74 | 84 | 92 |
 

<br>

Finally, this [table](https://drive.google.com/file/d/1NlkXbPSoMAATgeAqKidWvT6aU71wlpAf/view?usp=sharing) shows how many miRNAs are selected in at least one of the cross-validation. # column represents the number of cross-validation set of the corresponding miRNA from the table including avgerage(|Y_{tst} - Y_{pred}|).  


## Point3.b. prediction (Y) by tissue data (B) and cell line (X)

We trained the matrix B with the 240 cell lines. we call this B_{cell}.  
We applied B_{cell} to Y_{pred_ct} to see how the difference between the predicted expression and the true expression.  
- Y_{pred_ct} means predicted miRNA expression by B **c**ell line and X **t**issue.  
- Y_{tis} is miRNA expression of tissue data.  

![Image](/images/formula8.gif)    

The difference between Y_{pred_ct} and Y_{tis} are in [This table](https://drive.google.com/file/d/1KFGrGpgn_3YLOqT0Ar18TRAjZ2Lr8Bnf/view?usp=sharing).  

The three last column in the table,  
- **avg(diff)**: avgerage(|Y_{tis} - Y_{pred_ct}|)  
- **avg(distance)**: the average of a pair from 93 tissues. In this case, Y_{tis} has 93 tissues.  
Each cell has Y_{tis}, Y_{pred}, |Y_{tis} - Y_{pred_ct}|.  

**92 miRNAs** are meets avg(diff)<avg(distance).  

<br>
In addition, we trained the matrix B with the new tissue data; B_{tis} and then calculated |Y_{pred_tt} - Y_{tis}|, where Y_{pred_tt} the predicted expression and Y_{tis} is  true expression of miRNAs across the tissues.  

- distance_ct = |Y_{tis} - Y_{pred_ct}|
- distance_tt = |Y_{pred_tt} - Y_{tis}|
- Y_{pred_tt} means predicted miRNA expression by B **t**issue line and X **t**issue

So we compared distance_ct to distance_tt.  
[This shows](https://drive.google.com/file/d/1TOTZeGrePZCZqvFZKKrm37BIr8rNlM-9/view?usp=sharing) is the result.
Each cell has distance_ct, distance_tt, |distance_ct - distance_tt|.

**74 miRNAs** are satisfied with avg(distance_ct)<avg(distance_tt).  
**138 miRNAs** satisfy either avg(diff)< avg(distance) or avg(distance_ct)<avg(distance_tt).  

We also measured the absolute difference of the predicted expression and the true expression of a miRNA in a tissue compared with the absolute value of the true expression of this miRNA in this tissue.  
Y_{ratio} = |Y_{pred_ct} - Y_{trn}| / Y_{trn}  
This [table](https://drive.google.com/file/d/1BUZDMlkq6Z0L_K12XLHo7talvkaqZdem/view?usp=sharing) and a summary of the ratios.  

**107 miRNA** satisfy avg(diff_ratio) < avg(distance_ratio) from the above table.  

From the three sets; 92, 74, 107 miRNAs, I check how many overlapped each other.  
The below figure shows the result.  
![Image](/images/Figure9.gif)    


## Point3.c. the statistics of small distance
The avg(distance) and std(distance) are very similar. the correlation coefficient of two is 0.98 on the 3.b data, which means vert similar.  
The below figure shows the histogram of avg(distance). the y-axis is log-scale.    
most of them are in 200.  
![Image](/images/Figure7.png)    

The below figure is the histogram of avg(distance) < 200.  
![Image](/images/Figure8.png)  



